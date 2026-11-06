import datasketch
import copy
import logging
import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)

# Seed of the WeightedMinHashGenerator used by `build`
HASH_SEED = 1


def generate_treeweights( mastertree, taxaIndex ,  taxfilter= None, taxmask=None ):
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
    logger.info('making tree weights with %d taxa', len(taxaIndex))
    newtree = copy.deepcopy(mastertree)
    if taxmask:
        for node in newtree.traverse():
            if taxmask:
                if str(node.name) == str(taxmask):
                    logger.info('Masking tree to %s', taxmask)
                    newtree = node
                    break
        logger.info(
            'Masked tree contains %d taxa',
            len([n.name for n in newtree.traverse()]),
        )
    
    if taxfilter:
        for n in newtree.traverse():
            if n.name in taxfilter:            
                n.delete()
                
    for i in range(3):
        for n in newtree.traverse():
            weights[taxmax*i + taxaIndex[n.name] ] = 1
    return weights

def weighted_minhash(wmg, v):
    """
    Same result as `wmg.minhash(v)` (datasketch.WeightedMinHashGenerator),
    vectorized over the samples and restricted to the nonzero entries of v.

    datasketch loops over the samples and computes every one on the full
    dimension, with zeros masked as NaN; profiles are sparse, so this is
    much faster. The float32 arithmetic is the same, element by element, and
    the first minimum among nonzero entries is the one nanargmin picks.
    """
    v = np.asarray(v, dtype=np.float32)
    if v.shape != (wmg.dim,):
        raise ValueError("Input dimension mismatch, expecting %d" % wmg.dim)
    nonzero = np.flatnonzero(v)
    if nonzero.size == 0:
        raise ValueError("Input is all zeros")

    vlog = np.log(v[nonzero])
    rs = wmg.rs[:, nonzero]
    betas = wmg.betas[:, nonzero]
    t = np.floor((vlog / rs) + betas)
    ln_y = (t - betas) * rs
    ln_a = wmg.ln_cs[:, nonzero] - ln_y - rs
    k = np.argmin(ln_a, axis=1)

    hashvalues = np.zeros((wmg.sample_size, 2), dtype=int)
    hashvalues[:, 0] = nonzero[k]
    hashvalues[:, 1] = t[np.arange(wmg.sample_size), k]
    return datasketch.WeightedMinHash(wmg.seed, hashvalues)


def hash_tree(tp , taxaIndex , treeweights , wmg , lossonly = False , duplonly = False ):
    """
    Generate a weighted minhash and binary matrix row for a tree profile

    :param tp: a tree profile (pyham treemap or pyhamutils.TreeProfile)
    :param taxaIndex: dict mapping taxa to columns
    :param treeweights: a vector of weights for each tax levels
    :param wmg: Datasketch weighted minhash generator
    :return hog_matrix: a vector of weights for each tax level
    :return weighted_hash: a weighted minhash of a HOG

    """
    if not tp:
        return None, None

    taxaIndex_max = max(taxaIndex.values())+1
    hog_matrix_weighted = np.zeros(3*taxaIndex_max)
    hog_matrix_binary = np.zeros((1, 3*taxaIndex_max))

    presence, losses, dupl = [], [], []
    for n in tp.traverse():
        index = taxaIndex.get(n.name)
        if index is None:
            continue
        if n.nbr_genes > 0:
            presence.append(index)
        if n.lost:
            losses.append(index)
        if n.dupl:
            dupl.append(index)

    ### Athina comment:
    ### change things here to do permutation test!!!!
    for i, (event, indices) in enumerate((('presence', presence), ('loss', losses), ('dup', dupl))):
        if not indices:
            continue
        hogindex = np.asarray(indices) + i*taxaIndex_max

        if lossonly and event != 'loss':
            continue
        if duplonly and event != 'dup':
            continue
        hog_matrix_weighted[hogindex] = treeweights[hogindex, :].ravel()
        if not lossonly and not duplonly:
            hog_matrix_binary[:, hogindex] = 1

    if wmg.dim == hog_matrix_weighted.size and hog_matrix_weighted.sum() > 0:
        weighted_hash = weighted_minhash(wmg, hog_matrix_weighted)
        return  hog_matrix_binary , weighted_hash

    else:
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
        return {}

    hashes_subhogs = {}
    for key, tp in hog_tps['tree_dicts'].items():
        matrix, hashed = hash_tree(tp, taxaIndex, treeweights, wmg,
                                   lossonly=lossonly, duplonly=duplonly)
        if hashed is not None:
            hashes_subhogs[key] = (matrix, hashed)
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

def fam2hash_hdf5(fam, hdf5, dataset=None, nsamples=128, seed=HASH_SEED):
    """
    Read the stored hash values and return a weighted minhash
    :param fam: hog id
    :param hdf5: h5py object of the hashvalues
    :param dataset: which dataset to use when constructing the hash
    :param nsamples: number of permutations the index was built with
    :param seed: seed of the generator the index was built with
    :return: minhash1: the weighted hash of your HOG
    """
    if dataset is None:
        dataset = next(iter(hdf5.keys()))

    hash_values = np.asarray(hdf5[dataset][fam, :].reshape(nsamples,2 )).astype('int64')
    minhash1 = datasketch.WeightedMinHash(seed=seed, hashvalues=hash_values)
    return minhash1
