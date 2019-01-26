import ete3
import pandas as pd
from Bio import Entrez
import pyoma
import copy
import pickle

import os


def get_tree(oma=None, saveTree=True, overwrite = False):
    """
    Generates a working species tree for all OMA
    :param oma:  a pyoma db object
    :param saveTree: Bool for whether or not to save a mastertree newick file
    :return: tree_string: a newick string tree: an ete3 object
    """
    if not os.path.isfile('./mastertree.pkl') or overwrite == True:

        ncbi = ete3.NCBITaxa()

        genomes = pd.DataFrame(oma.root.Genome.read())["NCBITaxonId"].tolist()
        genomes = [ str(g) for g in genomes]
        tax = genomes + [ 131567, 2759, 2157, 45596 ]+[ taxrel[0] for taxrel in  list(oma.root.Taxonomy[:]) ]  + [  taxrel[1] for taxrel in list(oma.root.Taxonomy[:]) ]
        #add luca
        #tree_string = pyoma.browser.db.Taxonomy(oma.root.Taxonomy[:]).newick()
        #with open( './pyoma.nwk' , 'w') as nwkout:
        #    nwkout.write(tree_string)
        #print(tree_string)
        #tree_string = ete3.Tree( tree_string , format=1 )
        #ncbi.update_taxonomy_database()

        tax = set(tax)
        genomes = set(genomes)

        tax.remove(0)
        print(len(tax))

        tree = ete3.PhyloTree( name = '')
        tree.add_child(name ='131567')

        topo = ncbi.get_topology(tax , collapse_subspecies=False)
        tax = set([ str(taxid) for taxid in tax])


        tree.add_child(topo)

        orphans = list(genomes - set([x.name for x in tree.get_leaves()]))

        print('missing taxa:')
        print(len(orphans))

        Entrez.email = "clement.train@gmail.com"
        orphans_info1 = {}
        orphans_info2 = {}

        for x in orphans:
            search_handle = Entrez.efetch('taxonomy', id=str(x), retmode='xml')
            record = next(Entrez.parse(search_handle))
            print(record)
            orphans_info1[ record['ParentTaxId']] = x
            orphans_info2[x] = [x['TaxId'] for x in record['LineageEx']]

        for n in tree.traverse():
            if n.name in orphans_info1:
                n.add_sister(name = orphans_info1[n.name])
                print(n)


        orphans = set(genomes) - set([x.name for x in tree.get_leaves()])
        print(orphans)
        tree = add_orphans(orphans_info2, tree, genomes)
        orphans = set(genomes) - set([x.name for x in tree.get_leaves()])
        print(orphans)
        tree_string = tree.write(format=1)

        if saveTree == True:
            with open( './mastertree.nwk' , 'w') as nwkout:
                nwkout.write(tree_string)
            with open( './mastertree.pkl', 'wb') as pklout:
                pklout.write(pickle.dumps(tree))

    else:
        tree = pickle.loads(open('./mastertree.pkl', 'rb').read())
        tree_string = tree.write(format=1)
    return tree_string, tree



def generate_taxa_index(tree , taxfilter, taxmask):
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
