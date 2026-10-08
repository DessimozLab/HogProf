import xml.etree.ElementTree as ET
import logging
import os
import sys
import pyham

logger = logging.getLogger(__name__)


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

def _root_treemap(ham, root_hog, *, dataset_nodes, **kwargs):
    """Make a treemap for the rootHOG"""
    if dataset_nodes is not None and root_hog.genome.name not in dataset_nodes:
        return None

    return ham.create_tree_profile(hog=root_hog).treemap


def _check_limits(treemap, limit_events):
    """
    Check if |losses| AND |duplications| are less than limit_events
    """
    duplications = sum(node.dupl or 0 for node in treemap.traverse())
    losses = sum(node.lost or 0 for node in treemap.traverse())
    return duplications > limit_events or losses > limit_events


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


def _subhog_treemaps(ham, root_hog, *, family, limit_species,
                     limit_events, dataset_nodes, verbose):
    """Make treemaps for the levels mode"""
    hogs = root_hog.get_all_descendant_hogs()

    fallback_id = root_hog.hog_id if root_hog.hog_id is not None else family
    treemaps = {}

    # for every subhog
    for key, hog in zip(_subhog_keys(hogs, fallback_id), hogs):

        # check that the subHOG taxa was not filtered out
        if dataset_nodes is not None and hog.genome.name not in dataset_nodes:
            continue

        # check the minimal number of represented species
        if len(hog.get_all_descendant_genes_clustered_by_species()) < limit_species:
            continue

        treemap = ham.create_tree_profile(hog=hog).treemap

        # check the minimal number of events
        if _check_limits(treemap, limit_events):
            treemaps[key] = treemap

    if verbose:
        logger.debug("Family %s: retained %d of %d HOG profiles",
                     family, len(treemaps), len(hogs))
    return treemaps


# Profiling strategies supported in HogProf
_PROFILE_STRATEGIES = {
    # the default strategy -- 1 profile per rootHOG
    "root": _root_treemap,

    # Levels -- 1 profile per rootHOG + every subHOG
    "levels": _subhog_treemaps
}


def get_ham_treemap_from_row(row, tree_string,
                             *,
                             swap_ids=True,
                             orthoXML_as_string=True,
                             use_internal_name=True,
                             limit_species=10, limit_events=0,
                             dataset_nodes=None, verbose=False,
                             slicesubhogs=False):
    """Parse one root-HOG input and produce treemaps depending on the mode:

    - Root mode (default) returns a treemap or None
    - Levels mode (slicesubhogs=True) returns a dictionary,
    including the root when it passes the species and event thresholds.
    """
    family, orthoxml = row

    if not orthoxml:
        logger.warning("Family %s: empty OrthoXML provided", family)
        return {} if slicesubhogs else None

    if swap_ids:
        if not orthoXML_as_string:
            with open(orthoxml) as source:
                orthoxml = source.read()
        orthoxml = switch_name_ncbi_id(orthoxml)
        orthoXML_as_string = True

    ham = pyham.Ham(tree_string, orthoxml,
                    type_hog_file="orthoxml",
                    tree_format="newick_string",
                    use_internal_name=use_internal_name,
                    orthoXML_as_string=orthoXML_as_string)

    roots = ham.get_list_top_level_hogs()
    if len(roots) != 1:
        raise ValueError(f"Family {family}: expected one top-level HOG, found {len(roots)}")

    # pick the profiling strategy
    strategy = _PROFILE_STRATEGIES["levels" if slicesubhogs else "root"]

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
