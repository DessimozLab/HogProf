import numpy as np

#import sys
#sys.path.append('../')

from tables import *
from utils import goatools_utils
from goatools.semantic import TermCounts
from goatools import obo_parser
from goatools.associations import read_gaf
import pickle
import time
from utils import config_utils
import h5py
from functools import partial
import sys
sys.setrecursionlimit(500000)

import multiprocessing as mp

class Validation_semantic_similarity(object):

    def __init__(self, go_file, go_terms , gaf, omadb= None, tarfile_ortho = None , TermCountsFile = None):
        self.go_file = go_file

        if omadb :
            print('open oma db obj')
            from pyoma.browser import db
            h5_oma = open_file(omadb, mode="r")
            self.db_obj = db.Database(h5_oma)
            print('done')
        elif tarfile_ortho:
            #retrieve hog members from tarfile_ortho
             self.tar = tarfile.open(tarfile_ortho , "r:gz")
        else:
            raise Exception('please provide input dataset')

        #go_terms_hdf5 = h5py.File(go_terms, mode='r')
        #self.goterms2parents = go_terms_hdf5['goterms2parents']
        self.godf = pickle.loads( open(go_terms, 'rb').read())
        self.go_file = obo_parser.GODag(go_file)
        print('building gaf')
        self.gaf = goatools_utils.buildGAF( gaf)
        print('done')
        if TermCountsFile is None:
            self.termcounts = TermCounts(self.go_file, self.gaf)
        else:
            self.termcounts = pickle.loads(open(TermCountsFile ,'rb').read())
        #make a partial
        self.resniksimpreconf = partial( goatools_utils.resnik_sim_pandas , df = self.godf,  termcounts = self.termcounts)


    def semantic_similarity_score(self, hog_id_1, hog_id_2):
        """
        Runs semantic similarity analysis from 2 hog ids
        :param hog_id_1: first hog id
        :param hog_id_2: second hog id
        :return: semantic similarity score between the two hog ids
        """

        go_terms_1 = goatools_utils.get_go_terms_gaf(hog_id_1, self.db_obj , self.gaf)
        go_terms_2 = goatools_utils.get_go_terms_gaf(hog_id_2, self.db_obj , self.gaf)

        score = self._compute_score(go_terms_1, go_terms_2)
        return score


    def semantic_similarity_score_mp(self, hog_id_1, hog_id_2 , retq,  lock):
        """
        Runs semantic similarity analysis from 2 hog ids
        :param hog_id_1: first hog id
        :param hog_id_2: second hog id
        :return: semantic similarity score between the two hog ids
        """

        lock.acquire()
        go_terms_1 = goatools_utils.get_go_terms_gaf(hog_id_1, self.db_obj , self.gaf)
        go_terms_2 = goatools_utils.get_go_terms_gaf(hog_id_2, self.db_obj , self.gaf)
        lock.release()
        score = self._compute_score(go_terms_1, go_terms_2)
        retq.put( ((hog_id_1, hog_id_2) , score) )

        return score


    def _compute_genes_distance_cheap(self, go_terms_genes_1, go_terms_genes_2):
        """
        Computes matrix of distance between the genes of two hogs
        :param go_terms_genes_1: dictionary of genes
        :param go_terms_genes_2: dictionary of genes
        :return: matrix of distance between genes of hogs
        """
        if type(go_terms_genes_1) is dict and type(go_terms_genes_2) is dict and len(go_terms_genes_1)>0 and len(go_terms_genes_2):
            gos1 = list(go_terms_genes_1.values())
            gos2 = list( go_terms_genes_2.values())
            if len(gos2)>0 and len(gos1)>0:
                try:
                    setgo1 = set(gos1[0]).union( *gos1 )
                    setgo2 = set(gos2[0]).union( *gos2 )
                    keys=[]
                    for go1 in setgo1:
                        for go2 in setgo2:
                            keys.append(tuple(sorted((go1,go2))))
                    #infocontent of each term
                    res = [ self.resniksimpreconf(tup) for tup in keys]
                    res = dict( zip(keys,res))
                    genedist = np.zeros((len(go_terms_genes_1), len(go_terms_genes_2)))
                    normalizedgenedist = np.zeros((len(go_terms_genes_1), len(go_terms_genes_2)))
                    maxinf =   {go:self.resniksimpreconf((go,go)) for go in setgo1.union(setgo2) }
                    for i,gene1 in enumerate(go_terms_genes_1):
                        for j,gene2 in enumerate(go_terms_genes_2):
                            keyset =set([])
                            #generate all possible keys
                            [ keyset.add(tuple(sorted((go1,go2)))) for go1 in go_terms_genes_1[gene1] for go2 in go_terms_genes_2[gene2] ]
                            #get all go terms for these two genes
                            unique = set( go_terms_genes_1[gene1].union( go_terms_genes_2[gene2] ) )
                            if len(keyset)>0:
                                genedist[i,j] = np.amax( [ res[gopair]  for gopair in keyset  ] )
                                normalizedgenedist[i,j] = genedist[i,j] / np.amax([ maxinf[go] for go in unique ])
                            else:
                                genedist[i,j] = 0
                    return genedist, normalizedgenedist
                except:
                    return -1 , -1
        else:
            return -1 , -1

    def _compute_score(self, query_go_terms, result_go_terms):
        """
        Computes semantic similarity score between two hogs
        :param query_go_terms: dict of genes: list of go terms
        :param result_go_terms: dict of genes: list of go terms
        :return: semantic similarity score
        """

        dist_mat, normalized = self._compute_genes_distance_cheap(query_go_terms, result_go_terms)
        if type(dist_mat) is int and dist_mat ==-1:
            return -1 , -1

        score = self._mean_max_score_matrix(dist_mat)
        nscore = self._mean_max_score_matrix(normalized)
        #print(score)
        #only return nscore for now..

        return  nscore, score


    def _compute_genes_distance(self, go_terms_genes_1, go_terms_genes_2):
        """
        Computes matrix of distance between the genes of two hogs
        :param go_terms_genes_1: dictionary of genes
        :param go_terms_genes_2: dictionary of genes
        :return: matrix of distance between genes of hogs
        """
        # try:
        if type(go_terms_genes_1) is dict and type(go_terms_genes_2) is dict:
            keys_1 = go_terms_genes_1.keys()
            keys_2 = go_terms_genes_2.keys()

            gene_dist = np.zeros((len(keys_1), len(keys_2)))
            for k in range(len(keys_1)):
                for l in range(len(keys_2)):

                    go_terms_1 = list(go_terms_genes_1[list(keys_1)[k]])
                    go_terms_2 = list(go_terms_genes_2[list(keys_2)[l]])

                    if go_terms_1 and go_terms_2:
                        gene_dist[k, l] = self._compute_go_terms_score_per_gene(go_terms_1, go_terms_2)
            # except:
            #     gene_dist = -1
            return gene_dist
        else:
            return -1

    def _compute_go_terms_score_per_gene(self, go_terms_gene_1, go_terms_gene_2):
        """
        Computes the semantic similarity score between two genes
        :param go_terms_gene_1: list of go terms from one of the gene
        :param go_terms_gene_2: list of go terms from one of the gene
        :return: semantic similarity score between two genes
        """

        ss_dist = np.zeros((len(go_terms_gene_1), len(go_terms_gene_2)))

        for m in range(len(go_terms_gene_1)):
            for n in range(len(go_terms_gene_2)):
                #dist = goatools_utils.resnik_sim_hdf5(go_terms_gene_1[m], go_terms_gene_2[n], self.go_file, self.termcounts, self.goterms2parents)
                dist = goatools_utils.resnik_sim_pandas(go_terms_gene_1[m], go_terms_gene_2[n],self.godf,  self.termcounts)
                ss_dist[m, n] = dist
        gene_score = self._mean_max_score_matrix(ss_dist)

        return gene_score

    @staticmethod
    def _mean_max_score_matrix(matrix):
        """
        Computes the BMA of a matrix
        :param matrix: matrix
        :return: score: BMA of the matrix; returns -1 if matrix has 0 or 1 dimension
        """
        matrix_size = np.prod(matrix.shape)
        if not matrix_size:
            return -1
        score =  np.sum( matrix[np.where( matrix > 0 )] )  / matrix_size

        return score
