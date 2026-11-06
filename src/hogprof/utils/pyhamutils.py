import xml.etree.ElementTree as ET
import contextlib
import logging
import os
import sys
import types
import pyham
import pyham.ham
import pyham.taxonomy

logger = logging.getLogger(__name__)


class _CachedTaxonomies:
    """
    Per-process cache of pyham taxonomies, one per tree string.

    pyham.Ham parses the newick species tree for every family, and looks species
    up by a linear search over all tree nodes. For a ~5000 node tree this is
    most of the cost of small families. Here the tree is parsed once, nodes are
    indexed by name, and the per-family state pyham attaches to the nodes
    (the "genome" feature) is cleared before every reuse.

    A cached taxonomy must not be shared by two live Ham objects: families are
    processed one after another, and treemaps hold no references to the tree.
    """

    def __init__(self):
        self._cache = {}

    def get(self, tree_string, **kwargs):
        key = (tree_string, tuple(sorted(kwargs.items())))
        taxonomy = self._cache.get(key)
        if taxonomy is None:
            taxonomy = pyham.taxonomy.Taxonomy(tree_string, **kwargs)
            self._index(taxonomy)
            self._cache[key] = taxonomy
        else:
            self._reset(taxonomy)
        return taxonomy

    @staticmethod
    def _index(taxonomy):
        tree = taxonomy.tree
        by_name = {}
        for node in tree.traverse():
            by_name.setdefault(node.name, []).append(node)
        search_nodes = tree.search_nodes

        def fast_search_nodes(**conditions):
            if len(conditions) == 1 and "name" in conditions:
                return list(by_name.get(conditions["name"], ()))
            return search_nodes(**conditions)

        # pyham only calls tree.search_nodes(name=...)
        tree.search_nodes = fast_search_nodes

    @staticmethod
    def _reset(taxonomy):
        for node in taxonomy.leaves | taxonomy.internal_nodes:
            if "genome" in node.features:
                node.del_feature("genome")
        taxonomy.leaves = set()
        taxonomy.internal_nodes = set()


_taxonomies = _CachedTaxonomies()


@contextlib.contextmanager
def _use_cached_taxonomy():
    """Make pyham.Ham take its Taxonomy from the per-process cache"""
    original = pyham.ham.tax
    pyham.ham.tax = types.SimpleNamespace(
        Taxonomy=lambda tree_file, **kwargs: _taxonomies.get(tree_file, **kwargs))
    try:
        yield
    finally:
        pyham.ham.tax = original


class ProfileNode:
    """One level of a tree profile, with the fields of a pyham TreeProfile node"""
    __slots__ = ("name", "nbr_genes", "dupl", "lost", "retained",
                 "duplication", "nbr_events", "gain")

    def __init__(self, name, nbr_genes, dupl, lost, retained, duplication, nbr_events):
        self.name = name
        self.nbr_genes = nbr_genes
        self.dupl = dupl
        self.lost = lost
        self.retained = retained
        self.duplication = duplication
        self.nbr_events = nbr_events
        self.gain = None


class TreeProfile:
    """
    Lightweight replacement for pyham's TreeProfile.treemap: the levels below a
    HOG with their gene counts and events, without the tree structure.

    Supports the parts of the ete3 treemap HogProf uses: `.name` (the level of
    the HOG) and `.traverse()`.
    """
    __slots__ = ("name", "nodes")

    def __init__(self, name, nodes):
        self.name = name
        self.nodes = nodes

    def traverse(self):
        return iter(self.nodes)


def tree_profile(hog, genes_by_species=None):
    """
    Compute the tree profile of a HOG: the equivalent of
    `ham.create_tree_profile(hog=hog).treemap` (pyham.TreeProfile.computeTP_hog).

    Unlike pyham, works on the taxonomy nodes directly instead of a copy of the
    subtree made through newick. Besides being much faster, that keeps node
    names intact: the newick copy replaces characters such as '(' and ')' with
    '_', so species like 'Guillardia theta (strain CCMP2712)' did not match
    their genes and the taxa index.
    """
    root = hog.genome.taxon
    nodes = list(root.traverse())

    # map node with related hogs/genes
    level_groups = {node.name: [] for node in nodes}

    # add all of subhog to the related level
    for subhog in hog.get_all_descendant_hogs():
        if subhog.genome.name in level_groups:
            level_groups[subhog.genome.name].append(subhog)

    # add empty extant genomes
    for leaf in root.iter_leaves():
        level_groups[leaf.name] = []

    # add genes to related extant genome
    if genes_by_species is None:
        genes_by_species = hog.get_all_descendant_genes_clustered_by_species()
    for species, genes in genes_by_species.items():
        level_groups[species.name] = genes

    profile = []
    for lvl in nodes:
        if lvl is not root:
            cpt_dupl = 0
            cpt_ident = 0
            set_dup_parent = set()
            for h in level_groups[lvl.name]:
                if h.arose_by_duplication != False:
                    cpt_dupl += 1
                    set_dup_parent.add(h.parent)
                else:
                    cpt_ident += 1

            cpt_lost = 0
            for h_up in level_groups[lvl.up.name]:
                if len(h_up.children) != 0:
                    found = False
                    for child in h_up.children:
                        if child.genome.taxon.name == lvl.name:
                            found = True
                    if found is False:
                        cpt_lost += 1
                else:
                    cpt_lost += 1

            cpt_duplication = max(0, cpt_dupl - len(set_dup_parent))
            nbr_ev = cpt_lost + cpt_duplication
        else:
            cpt_dupl = cpt_ident = cpt_lost = cpt_duplication = nbr_ev = None

        profile.append(ProfileNode(
            lvl.name, len(level_groups[lvl.name]), cpt_dupl, cpt_lost,
            cpt_ident, cpt_duplication, nbr_ev))

    return TreeProfile(root.name, profile)


def get_orthoxml_oma(fam, db_obj):
    orthoxml = db_obj.get_orthoxml(fam).decode()    
    return orthoxml

def get_orthoxml_tar(fam, tar):
    f = tar.extractfile(fam)
    if f is not None:
        return f.read()
    else:
        raise Exception( member + ' : not found in tarfile ')
    return orthoxml

def get_species_from_orthoxml(orthoxml):
    NCBI_taxid2name = {}
    root = ET.fromstring(orthoxml)
    for child in root:
        if 'species' in child.tag:
            NCBI_taxid2name[child.attrib['NCBITaxId']] = child.attrib['name']
    return NCBI_taxid2name

def switch_name_ncbi_id(orthoxml , mapdict = None  ):
    #swap ncbi taxid for species name to avoid ambiguity
    #mapdict should be a mapping from species name to taxid if the info isnt in the orthoxmls
    root = ET.fromstring(orthoxml)
    for child in root:
        if 'species' in child.tag:
            child.attrib['name'] = child.attrib['NCBITaxId']
        elif mapdict:
            child.attrib['name'] = mapdict[child.attrib['name']]
        
    orthoxml = ET.tostring(root, encoding='unicode', method='xml')
    return orthoxml

def reformat_treenames( tree , mapdict = None  ):   
    #tree is an ete3 tree instance
    #replace ( ) - / . and spaces with underscores
    #iterate over all nodes
    for node in tree.traverse():
        if mapdict:
            node.name = mapdict[node.name]
        else:
            node.name = node.name.replace('(', '').replace(')', '').replace('-', '_').replace('/', '_')
    return tree

def reformat_names_orthoxml(orthoxml , mapdict = None  ):
    #replace ( ) - / . and spaces with underscores
    root = ET.fromstring(orthoxml)
    for child in root:
        if 'species' in child.tag:
            child.attrib['name'] = child.attrib['name'].replace('(', '').replace(')', '').replace('-', '_').replace('/', '_')
        elif mapdict:
            child.attrib['name'] = mapdict[child.attrib['name']]
    orthoxml = ET.tostring(root, encoding='unicode', method='xml')
    return orthoxml

def create_nodemapping(tree):
    #create a mapping from node name to node
    nodemapping = {}
    for i,node in enumerate(tree.traverse()):
        nodemapping[node.name] = str(i)
    
    #assert number of nodes is equal to number of unique names
    assert len([n for n in tree.traverse()]) == len(set(nodemapping.values()))

    return nodemapping

def tree2numerical(tree):
    mapper = create_nodemapping(tree)
    for i,node in enumerate(tree.traverse()):
        node.name = mapper[node.name]
    return tree , mapper

def orthoxml2numerical(orthoxml , mapper):
    #print(orthoxml)
    try:
        root = ET.fromstring(orthoxml)
    except Exception as e:
        with open( orthoxml , 'r') as f:
            orthoxml = f.read()
        root = ET.fromstring(orthoxml)
    for child in root:
        if 'name' in child.attrib:
            child.attrib['name'] = mapper[child.attrib['name']]
    orthoxml = ET.tostring(root, encoding='unicode', method='xml')
    return orthoxml 

def _root_treemap(ham, root_hog, *, dataset_nodes, limit_events, **kwargs):
    """Make a treemap for the rootHOG"""
    if dataset_nodes is not None and root_hog.genome.name not in dataset_nodes:
        return None

    treemap = tree_profile(root_hog)
    return treemap if _has_min_events(treemap, limit_events) else None


def count_events(treemap):
    """Total number of losses and duplications in a tree profile"""
    return sum((node.dupl or 0) + (node.lost or 0) for node in treemap.traverse())


def _has_min_events(treemap, limit_events):
    """
    Check that the profile has at least `limit_events` events
    (losses + duplications). `limit_events=0` keeps every profile.
    """
    return limit_events <= 0 or count_events(treemap) >= limit_events


def _subhog_keys(hogs, family):
    """Generate keys for subhogs in family"""
    bases = [
        f"{hog.genome.name}_{hog.hog_id if hog.hog_id is not None else family}"
        for hog in hogs
    ]
    reserved = set(bases)
    seen = set()
    suffixes = {}
    for base in bases:
        key = base
        suffix = suffixes.get(base, 2)
        while key in seen or (key != base and key in reserved):
            key = f"{base}__{suffix}"
            suffix += 1
        suffixes[base] = suffix
        seen.add(key)
        yield key


def _iter_subhog_treemaps(ham, root_hog, *, family, limit_species,
                          limit_events, dataset_nodes, verbose):
    """
    Yield (key, treemap) for the levels mode, one at a time.

    A treemap is a copy of the species subtree below the subHOG, so a large
    family (tens of thousands of subHOGs) needs many GB if all are kept at once.
    Consumers should hash each treemap and drop it.
    """
    hogs = root_hog.get_all_descendant_hogs()

    fallback_id = root_hog.hog_id if root_hog.hog_id is not None else family
    retained = 0

    # for every subhog
    for key, hog in zip(_subhog_keys(hogs, fallback_id), hogs):

        # check that the subHOG taxa was not filtered out
        if dataset_nodes is not None and hog.genome.name not in dataset_nodes:
            continue

        # check the minimal number of represented species
        genes_by_species = hog.get_all_descendant_genes_clustered_by_species()
        if len(genes_by_species) < limit_species:
            continue

        treemap = tree_profile(hog, genes_by_species)

        # check the minimal number of events
        if _has_min_events(treemap, limit_events):
            retained += 1
            yield key, treemap

    if verbose:
        logger.debug("Family %s: retained %d of %d HOG profiles",
                     family, retained, len(hogs))


def _subhog_treemaps(ham, root_hog, **kwargs):
    """Make treemaps for the levels mode"""
    return dict(_iter_subhog_treemaps(ham, root_hog, **kwargs))


# Profiling strategies supported in HogProf
_PROFILE_STRATEGIES = {
    # the default strategy -- 1 profile per rootHOG
    "root": _root_treemap,

    # Levels -- 1 profile per rootHOG + every subHOG
    "levels": _subhog_treemaps,

    # Levels, lazily: yields (key, treemap) pairs
    "levels-lazy": _iter_subhog_treemaps,
}


def get_ham_treemap_from_row(row, tree_string,
                             *,
                             swap_ids=True,
                             orthoXML_as_string=True,
                             use_internal_name=True,
                             limit_species=10, limit_events=0,
                             dataset_nodes=None, verbose=False,
                             slicesubhogs=False, lazy=False):
    """Parse one root-HOG input and produce treemaps depending on the mode:

    - Root mode (default) returns a treemap or None
    - Levels mode (slicesubhogs=True) returns a dictionary,
    including the root when it passes the species and event thresholds.
    With lazy=True it returns an iterator of (key, treemap) pairs instead.
    """
    family, orthoxml = row

    if not orthoxml:
        logger.warning("Family %s: empty OrthoXML provided", family)
        return ({} if not lazy else iter(())) if slicesubhogs else None

    if swap_ids:
        if not orthoXML_as_string:
            with open(orthoxml) as source:
                orthoxml = source.read()
        orthoxml = switch_name_ncbi_id(orthoxml)
        orthoXML_as_string = True

    with _use_cached_taxonomy():
        ham = pyham.Ham(tree_string, orthoxml,
                        type_hog_file="orthoxml",
                        tree_format="newick_string",
                        use_internal_name=use_internal_name,
                        orthoXML_as_string=orthoXML_as_string)

    roots = ham.get_list_top_level_hogs()
    if len(roots) != 1:
        raise ValueError(f"Family {family}: expected one top-level HOG, found {len(roots)}")

    # pick the profiling strategy
    if slicesubhogs:
        strategy = _PROFILE_STRATEGIES["levels-lazy" if lazy else "levels"]
    else:
        strategy = _PROFILE_STRATEGIES["root"]

    return strategy(ham, roots[0],
                    family=family,
                    limit_species=limit_species, limit_events=limit_events,
                    dataset_nodes=None if dataset_nodes is None else set(dataset_nodes),
                    verbose=verbose)


def get_subhog_ham_treemaps_from_row(row, tree_string, **kwargs):
    """Compatibility entry point for levels-mode callers."""
    return get_ham_treemap_from_row(row, tree_string, slicesubhogs=True, **kwargs)


def add_library_path(library_path):
    """Add the directory containing profiler.py to Python's sys.path"""
    if not os.path.exists(library_path):
        raise FileNotFoundError(f"Path does not exist: {library_path}")
    profiler_dir = os.path.dirname(os.path.abspath(library_path))
    if profiler_dir not in sys.path:
        sys.path.append(profiler_dir)


def yield_families(h5file, start_fam):
    """
    Given a h5file containing OMA server, returns an iterator over the families
    (not sure if still in use)
    :param h5file: omafile
    :param start_fam: fam to start on
    :return: fam number
    """
    for row in h5file.root.OrthoXML.Index:
        if row[0] > start_fam:
            yield row[0]
def get_one_family(i, h5file):
    '''
    get one family from database
    Args:
        i : family number
        h5file : OMA server file
    Return :
        family
    Not sure if still in use
    '''
    return h5file.root.OrthoXML.Index[i][0]
