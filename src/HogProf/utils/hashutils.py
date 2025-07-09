

import datasketch
import itertools
import ete3
import copy
import math
import numpy as np
import pandas as pd


def generate_treeweights( mastertree, taxaIndex ,  taxfilter, taxmask ):
    #weighing function for tax level, masking levels etc. sets all weights to 1 if they are in taxmask or filter
    #custom weights can also be used here
    """
    Generate the weights of each taxonomic level to be applied during the
    constructin of weighted minhashes
    :param mastertree: full corrected ncbi taxonomy
    :param taxaIndex: dict mapping taxa to columns
    :param taxfilter: list of branches to delete
    :param taxmask: if this is not NONE taxmask, the DB is constructed with this subtree
    :return: weights: a vector of weights for each tax level
    """
    #get max of taxa index
    taxmax = max(taxaIndex.values())+1
    weights = np.zeros((3*taxmax,1))
    print('making tree weights w n taxa = :',len(taxaIndex))
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
    for i in range(3):
        for n in newtree.traverse():
            weights[taxmax*i + taxaIndex[n.name] ] = 1
    return weights

def hash_tree(tp , taxaIndex , treeweights , wmg , lossonly = False , duplonly = False ):
    """
    Generate a weighted minhash and binary matrix row for a tree profile

    :param tp: a pyham tree profile
    :param taxaIndex: dict mapping taxa to columns
    :param treeweights: a vector of weights for each tax levels
    :param wmg: Datasketch weighted minhash generator
    :return hog_matrix: a vector of weights for each tax level
    :return weighted_hash: a weighted minhash of a HOG

    """
    if not tp:
        return None, None

    taxaIndex_max = max(taxaIndex.values())+1
    hog_matrix_weighted = np.zeros((1, 3*taxaIndex_max))
    hog_matrix_binary = np.zeros((1, 3*taxaIndex_max))

    losses = [ taxaIndex[n.name]  for n in tp.traverse() if n.lost and n.name in taxaIndex  ]
    dupl = [ taxaIndex[n.name]  for n in tp.traverse() if n.dupl  and n.name in taxaIndex  ]
    presence = [ taxaIndex[n.name]  for n in tp.traverse() if n.nbr_genes > 0  and n.name in taxaIndex ]
    indices = dict(zip (['presence', 'loss', 'dup'],[presence,losses,dupl] ) )
    ### Athina comment:
    ### change things here to do permutation test!!!!
    for i,event in enumerate(indices):
        if len(indices[event])>0:
            try:
                hogindex = np.asarray(indices[event])+i*taxaIndex_max
                
                hog_matrix_weighted[:,hogindex] = treeweights[hogindex , : ].ravel()
                
                if lossonly == True and event != 'loss':
                    hog_matrix_weighted[:,hogindex] = 0
                if duplonly == True and event != 'dup':
                    hog_matrix_weighted[:,hogindex] = 0
                if lossonly == False and duplonly == False:
                    hog_matrix_binary[:,hogindex] = 1
            except:
                print( 'error in hash_tree')
                print( 'event', event)
                print( 'indices', indices[event])
                print( 'hogindex', hogindex)

    input_vec = list(hog_matrix_weighted.flatten())

    if wmg.dim == len(input_vec):
        weighted_hash = wmg.minhash(input_vec)
        return  hog_matrix_binary , weighted_hash

    else:
        print('error in hash_tree')
        print('wmg.dim', wmg.dim)
        print('len(input_vec)', len(input_vec))
        print( input_vec)
        return None, None




def hash_trees_subhogs(hog_tps , taxaIndex , treeweights , wmg , lossonly = False , duplonly = False ):
    """
    Generate a weighted minhash and binary matrix row for multiple tree profiles within a HOG

    :param tp: a pyham tree profile
    :param taxaIndex: dict mapping taxa to columns
    :param treeweights: a vector of weights for each tax levels
    :param wmg: Datasketch weighted minhash generator
    :return hog_matrix: a vector of weights for each tax level
    :return weighted_hash: a weighted minhash of a HOG

    """
    #print(hog_tps.keys())
    #print(hog_tps[['tree_dicts']].to_dict())
    if hog_tps is None or 'tree_dicts' not in hog_tps or hog_tps['tree_dicts'] is None:
        return None

    hashes_subhogs = {key:hash_tree(tp , taxaIndex , treeweights , wmg) for key,tp  in hog_tps['tree_dicts'].items()}
    
    return hashes_subhogs


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

def row2hash(row , taxaIndex , treeweights , wmg , lossonly = False , duplonly = False):
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
    #print(f'Processing row in row2hash: fam={fam}, treemap={treemap}')
    hog_matrix,weighted_hash = hash_tree(treemap , taxaIndex , treeweights , wmg , lossonly = lossonly , duplonly = duplonly)
    #print(f'Generated hog_matrix: {hog_matrix}, weighted_hash: {weighted_hash}')
    
    return  pd.Series([weighted_hash,hog_matrix], index=['hash','rows'])

def fam2hash_hdf5(fam,  hdf5, dataset = None, nsamples = 128, fam2orthoxmlpath = None):
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
    ### Athina comment: this here below may need fixing for levels !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ### because fam works as an index only when slicesubhogs=False
    ### otherwise we need to get to the id again
    #print(hdf5[dataset][:])
    #print(fam)
    if fam2orthoxmlpath is None:
        hashvalues = np.asarray(hdf5[dataset][fam, :].reshape(nsamples,2 ))
    else:
        id2famsubhog_df = pd.read_csv(fam2orthoxmlpath, index_col=0)
        print(id2famsubhog_df.head())
        ### Athina note: this part takes a long time! (or not) - may not be working. {1383:[1383],1413:[1413],...}
        fam_dict = id2famsubhog_df.groupby('fam').apply(lambda x: x.index.tolist()).to_dict()
        #print(fam_dict)
        #subhog_dict = id2famsubhog_df.set_index('subhog_id').to_dict(orient='index')
        if isinstance(fam, int):
            indices = fam_dict[fam]
        elif isinstance(fam, list):
            indices = fam
        elif isinstance(fam, str):
            fam_parts = fam.split('_')
            fam_int = int(fam_parts[0])
            subhog_str = '_'.join(fam_parts[1:])
            indices = id2famsubhog_df[(id2famsubhog_df['fam'] == fam_int) & (id2famsubhog_df['subhog_id'] == subhog_str)].index.tolist() 
        print('recursive call to fam2hash_hdf5')
        return [fam2hash_hdf5(i, hdf5, dataset, nsamples) for i in indices]
    #print(hashvalues)
    hashvalues = hashvalues.astype('int64')
    minhash1 = datasketch.WeightedMinHash( seed = 1, hashvalues=hashvalues)
    return minhash1

