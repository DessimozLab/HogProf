

import datasketch
import itertools
import ete3
import copy
import math
import numpy as np
import pandas as pd


def generate_treeweights( mastertree, taxaIndex , taxfilter, taxmask ):
    #weighing function for tax level, masking levels etc. sets all weights to 1
    """
    Generate the weights of each taxonomic level to be applied during the
    constructin of weighted minhashes
    :param mastertree: full corrected ncbi taxonomy
    :param taxaIndex: dict mapping taxa to columns
    :param taxfilter: list of branches to delete
    :param taxmask: if this is not NONE taxmask, the DB is constructed with this subtree
    :return: weights: a vector of weights for each tax level
    """

    weights = { type: np.zeros((len(taxaIndex),1)) for type in ['presence', 'loss', 'dup']}
    print(len(taxaIndex))
    newtree = mastertree
    for event in weights:
        for n in newtree.traverse():
            if taxmask:
                if str(n.name) == str(taxmask):
                    newtree = n
                    break
            if taxfilter:
                if n.name in taxfilter:
                    n.delete()
    for event in weights:
        for n in newtree.traverse():
            weights[event][taxaIndex[n.name]] = 1
    return weights

def hash_tree(tp , taxaIndex , treeweights , wmg):
    """
    Generate a weighted minhash and binary matrix row for a tree profile

    :param tp: a pyham tree profile
    :param taxaIndex: dict mapping taxa to columns
    :param treeweights: a vector of weights for each tax levels
    :param wmg: Datasketch weighted minhash generator
    :return hog_matrix: a vector of weights for each tax level
    :return weighted_hash: a weighted minhash of a HOG

    """

    losses = [ taxaIndex[n.name]  for n in tp.traverse() if n.lost and n.name in taxaIndex  ]
    dupl = [ taxaIndex[n.name]  for n in tp.traverse() if n.dupl  and n.name in taxaIndex  ]
    presence = [ taxaIndex[n.name]  for n in tp.traverse() if n.nbr_genes > 0  and n.name in taxaIndex ]
    indices = dict(zip (['presence', 'loss', 'dup'],[presence,losses,dupl] ) )
    hog_matrix_weighted = np.zeros((1, 3*len(taxaIndex)))
    hog_matrix_binary = np.zeros((1, 3*len(taxaIndex)))
    for i,event in enumerate(indices):
        if len(indices[event])>0:
            taxindex = np.asarray(indices[event])
            hogindex = np.asarray(indices[event])+i*len(taxaIndex)
            hog_matrix_weighted[:,hogindex] = treeweights[hogindex].ravel()
            hog_matrix_binary[:,hogindex] = 1
    weighted_hash = wmg.minhash(list(hog_matrix_weighted.flatten()))

    return  hog_matrix_binary , weighted_hash

def tree2str_DCA(tp , taxaIndex ):
    """
    Generate a string where each column is a tax level
    each letter code corresponds to an event type
    each row is a protein family. for use with DCA pipelines

    :param tp: a pyham tree profile
    :param taxaIndex: dict mapping taxa to columns
    :return dcaMat: a weighted minhash of a HOG
    """
    #convert a tree profile to a weighted minhash


    losses = [ taxaIndex[n.name]  for n in tp.traverse() if n.lost and n.name in taxaIndex  ]
    dupl = [ taxaIndex[n.name]  for n in tp.traverse() if n.dupl  and n.name in taxaIndex  ]
    presence = [ taxaIndex[n.name]  for n in tp.traverse() if n.nbr_genes > 0  and n.name in taxaIndex ]


    Ds = list( set(dupl).intersection(set(presence)))
    Ps=  list(set(presence).difference(set(dupl)))
    Ls=  list(set(losses))
    charar = np.chararray(len(taxaIndex) )
    #set to absent

    charar.fill( 'A')
    charar[Ds] = 'D'
    charar[Ls] = 'L'
    charar[Ps] = 'P'

    return charar

def row2hash(row , taxaIndex , treeweights , wmg):
    """
    turn a dataframe row with an orthoxml file to hash and matrix row
    :param row: lsh builder dataframe row
    :param taxaIndex: dict mapping taxa to columnsfam
    :param treeweights: a vector of weights for each tax levels
    :param wmg: Datasketch weighted minhash generator
    :return: hog_matrix: a vector of weights for each tax level
    :return: weighted_hash: a weighted minhash of a HOG
    """
    #convert a dataframe row to a weighted minhash
    fam, treemap = row.tolist()
    hog_matrix,weighted_hash = hash_tree(treemap , taxaIndex , treeweights , wmg)
    return [weighted_hash,hog_matrix]


def fam2hash_hdf5(fam,  hdf5, dataset = None, nsamples = 128  ):
    #read the stored hash values and return a weighted minhash
    """
    Read the stored hash values and return a weighted minhash
    :param fam: hog id
    :param hdf5: h5py object of the hashvalues
    :param dataset: which dataset to use when constructing the hash
    :return: minhash1: the weighted hash of your HOG
    """
    if dataset is None:
        dataset = list(hdf5.keys())[0]
    hashvalues = np.asarray(hdf5[dataset][fam, :].reshape(nsamples,2 ))
    hashvalues = hashvalues.astype('int64')
    minhash1 = datasketch.WeightedMinHash( seed = 1, hashvalues=hashvalues)
    return minhash1

def hogid2fam(hog_id):
    """
    For use with OMA HOGs
    Get fam given hog id
    :param hog_id: hog id
    :return: fam

    """

    if not hog_id:
        return hog_id
    if type(hog_id) is int:
        return hog_id

    if ':' in hog_id:
        hog_id = hog_id.split(':')[1]
        if '.' in hog_id:
            hog_id = hog_id.split('.')[0]
        hog_id = hog_id.replace("'",'')
        fam = int(hog_id)
    else:
        fam = int(hog_id)
    return fam


def fam2hogid(fam_id):
    """
    For use with OMA HOGs
    Get hog id given fam
    :param fam_id: fam
    :return: hog id
    """
    hog_id = "HOG:" + (7-len(str(fam_id))) * '0' + str(fam_id)
    return hog_id
