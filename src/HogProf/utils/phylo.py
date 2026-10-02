"""
Wrapper module isolating the ete3 backend.
Provides functions to read and load newick trees
and owns the details about the underlying format
that HogProf expects.

This module is expected to be used instead
directly importing ete3 and any of its functions.
"""

import copy
import pickle
import logging
from pathlib import Path
from typing import Union, Iterable
import ete3
from Bio import Entrez

logger = logging.getLogger(__name__)


# contains the key-value args passed to ete3 every time
# we need to parse the tree from newick
default_kwargs = {
    "format" : 1,
    "quoted_node_names" : True
}

def from_file(filename: Union[str, Path], **kwargs) -> ete3.Tree:
    logger.debug(f"Reading file {filename}")
    ete3_kwargs = default_kwargs.copy()
    ete3_kwargs.update(kwargs)
    file_str = Path(filename).absolute().as_posix()
    return ete3.Tree(file_str, **ete3_kwargs)

def from_string(string: str, **kwargs) -> ete3.Tree:
    ete3_kwargs = default_kwargs.copy()
    ete3_kwargs.update(kwargs)
    return ete3.Tree(string, **ete3_kwargs)

def to_file(tree: ete3.Tree, filename: Union[str, Path], **kwargs) -> None:
    logger.debug(f"Saving to file {filename}")
    ete3_kwargs = default_kwargs.copy()
    ete3_kwargs.update(kwargs)
    with open(filename, 'w') as outfile:
        outfile.write(tree.write(**ete3_kwargs))

def to_string(tree: ete3.Tree, **kwargs) -> str:
    ete3_kwargs = default_kwargs.copy()
    ete3_kwargs.update(kwargs)
    return tree.write(**ete3_kwargs)


def add_orphans(orphan_info, tree, genome_ids_list, verbose=False):
    """
    Fix the NCBI taxonomy by adding missing species.
    :param: orphan_info: a dictionary containing info from the NCBI on the missing taxa
    :param: tree : an ete3 tree missing some species
    :param: genome_ids_list: the comlete set of taxids that should be on the tree
    :verbose: Bool print debugging stuff
    :return: tree: a species tree with the orphan genomes added
    """
    first = True


    newdict = {}

    leaves = set([leaf.name for leaf in tree.get_leaves()])

    orphans = set(genome_ids_list) - leaves
    oldkeys = set(list(newdict.keys()))

    keys = set()
    i = 0
    print(i)

    while first or ( len(orphans) > 0  and keys != oldkeys ) :
        first = False
        oldkeys = keys
        leaves = set([leaf.name for leaf in tree.get_leaves()])
        orphans = set(genome_ids_list) - leaves
        print(len(orphans))
        for orphan in orphans:
            if str(orphan_info[orphan][-1]) in newdict:
                newdict[str(orphan_info[orphan][-1])].append(orphan)
            else:
                newdict[str(orphan_info[orphan][-1])] = [orphan]
        keys = set(list(newdict.keys()))
        for n in tree.traverse():
            if n.name in newdict and n.name not in leaves:
                for orph in newdict[n.name]:
                    n.add_sister(name=orph)
                del newdict[n.name]

        for orphan in orphans:
            if len(orphan_info[orphan]) > 1:
                orphan_info[orphan].pop()

        newdict = {}
    nodes = {}
    print(orphans)
    #clean up duplicates
    for n in tree.traverse():
        if n.name not in nodes:
            nodes[ n.name] =1
        else:
            nodes[ n.name] +=1

    for n in tree.traverse():
        if nodes[ n.name] >1:
            if n.is_leaf()== False:
                n.delete()
                nodes[ n.name]-= 1


    return tree


def get_tree(genomes, outdir=None):
    """
    Generates a taxonomic tree using the ncbi taxonomy and
    :param oma:  a pyoma db object
    :param saveTree: Bool for whether or not to save a mastertree newick file
    :return: tree_string: a newick string tree: an ete3 object

    """
    ncbi = ete3.NCBITaxa()
    genomes = set(genomes)
    tree = ete3.PhyloTree(name='-1')
    topo = ncbi.get_topology(genomes, collapse_subspecies=False)
    tree.add_child(topo)
    orphans = list(genomes - set([x.name for x in tree.get_leaves()]))
    print('missing taxa:')
    print(len(orphans))

    orphans_info1 = {}
    orphans_info2 = {}
    for x in orphans:
        search_handle = Entrez.efetch('taxonomy', id=str(x), retmode='xml')
        record = next(Entrez.parse(search_handle))
        print(record)
        orphans_info1[record['ParentTaxId']] = x
        orphans_info2[x] = [x['TaxId'] for x in record['LineageEx']]
    for n in tree.traverse():
        if n.name in orphans_info1:
            n.add_sister(name=orphans_info1[n.name])
            print(n)
    orphans = set(genomes) - set([x.name for x in tree.get_leaves()])
    tree = add_orphans(orphans_info2, tree, genomes)
    orphans = set(genomes) - set([x.name for x in tree.get_leaves()])
    tree_string = tree.write(format=1)

    with open(outdir + 'master_tree.nwk', 'w') as nwkout:
        nwkout.write(tree_string)
    with open(outdir + '_master_tree.pkl', 'wb') as pklout:
        pklout.write(pickle.dumps(tree))

    return tree_string, tree


def generate_taxa_index(tree , taxfilter= None, taxmask=None):
    """
    Generates an index for the global taxonomic tree for all OMA
    :param tree: ete3 tree
    :return: taxaIndex: dictionary key: node name (species name); value: index
        taxaIndexReverse: dictionary key: index: value: species name
    """
    newtree = copy.deepcopy(tree)
    for n in newtree.traverse():
        if taxmask:
            if str(n.name) == str(taxmask):
                newtree = n
                break
        if taxfilter:
            if n.name in taxfilter:
                #set weight for descendants of n to 0
                n.delete()
    taxa_index = {}
    taxa_index_reverse = {}
    for i, n in enumerate(tree.traverse()):
        taxa_index_reverse[i] = n.name
        taxa_index[n.name] = i-1

    return taxa_index, taxa_index_reverse


class TreeValidator:
    """Validate and repair a species tree before starting worker processes."""

    def __init__(self, filename: Path, species_names: Iterable[str] = ()):
        self.filename = Path(filename)
        self.species_names = {str(name) for name in species_names}

    def run(self):
        self._validate_format()
        tree = from_file(self.filename)

        self._fix_internal_species(tree)

        leaf_names_before = tree.get_leaf_names()
        new_tree = self._promote_single_leaves(tree.copy())

        # repaired tree check: kept every leaf
        assert sorted(new_tree.get_leaf_names()) == sorted(leaf_names_before)
        # check: new tree is connected
        assert all(child.up is node for node in new_tree.traverse() for child in node.children)

        return new_tree

    def _validate_format(self):
        # validate name
        valid = self.filename.suffix.lower() in [".nwk", ".newick"]
        if not valid:
            raise ValueError("Input tree must be in the newick format")

    def _check_missing_duplicated_species(self, tree, nodes_by_name):
        if not self.species_names:
            return

        # check if species list has the names not existing in the tree
        missing = self.species_names - nodes_by_name.keys()
        duplicated = {
            name for name, nodes in nodes_by_name.items() if len(nodes) != 1
        }
        if missing or duplicated:
            problems = []
            if missing:
                problems.append("missing species: " + ", ".join(sorted(missing)))
            if duplicated:
                problems.append(
                    "species mapped to multiple tree nodes: "
                    + ", ".join(sorted(duplicated))
                )
            raise ValueError("Invalid species tree (" + "; ".join(problems) + ")")


    def _fix_internal_species(self, tree):
        """
        Pyham raises the following error
        `TypeError: species name 'XXX' maps to an ancestral name, not a leaf of the taxonomy`
        if an internal node is marked as species. Repair the input tree for these cases
        """
        if not self.species_names:
            return

        nodes_by_name = {}
        for node in tree.traverse():
            if node.name in self.species_names:
                nodes_by_name.setdefault(node.name, []).append(node)

        self._check_missing_duplicated_species(tree, nodes_by_name)

        repaired = []
        # Postorder also handles the unlikely case of nested species nodes
        for node in tree.traverse("postorder"):
            if node.name not in self.species_names or node.is_leaf():
                continue
            if node.is_root():
                raise ValueError(
                    f"Species {node.name!r} maps to the root of a non-trivial tree"
                )

            parent = node.up
            for child in list(node.get_children()):
                child.detach()
                parent.add_child(child)
            repaired.append(node.name)

        if repaired:
            logger.warning(
                "Converted %d internal species nodes to leaves: %s",
                len(repaired),
                ", ".join(sorted(repaired)),
            )

        invalid = [
            name
            for name, nodes in nodes_by_name.items()
            if len(nodes) != 1 or not nodes[0].is_leaf()
        ]
        if invalid:
            raise ValueError(
                "Species did not resolve to leaves after tree repair: "
                + ", ".join(sorted(invalid))
            )

    @staticmethod
    def _promote_single_leaves(tree: ete3.Tree) -> ete3.Tree:
        """
        If a node has a single descendant that's a leaf,
        makes it a sister node
        """
        new_tree = tree

        count = 0
        for node in new_tree.traverse("postorder"):
            children = node.get_children()
            to_delete = not node.is_root() and len(children) == 1 and children[0].is_leaf()
            if to_delete:
                logger.debug("Reparenting node %s", children[0].name)
                count += 1
                node.delete()

        if count > 0:
            logger.info("Repaired %d internal single-child nodes", count)
        return new_tree

