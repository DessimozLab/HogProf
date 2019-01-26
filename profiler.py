
import _pickle as pickle
import pandas as pd
import h5py
import itertools
import ujson as json
import random
from scipy.sparse import csr_matrix
from tables import *
from pyoma.browser import db
import numpy as np
import random
np.random.seed(0)
random.seed(0)
import ete3
from datasketch import WeightedMinHashGenerator
from validation import validation_semantic_similarity
from utils import hashutils,  config_utils , pyhamutils , files_utils
from time import time
import multiprocessing as mp
import functools
import numpy as np
import time
import sys
import gc
import logging


class Profiler:


    def __init__(self,lshforestpath = None, hashes_h5=None, mat_path= None, oma = False , nsamples = 256):
        #use the lsh forest or the lsh

        """
        A profiler object allows the user to query the LSH with HOGs and get a list of result HOGs back

        """
        print('loading lsh')
        with open(lshforestpath, 'rb') as lshpickle:
            self.lshobj = pickle.loads(lshpickle.read())
            print('indexing lsh')
            self.lshobj.index()
        self.hashes_h5 = h5py.File(hashes_h5, mode='r')
        self.nsamples = nsamples
        print('DONE')

        if mat_path:
            ## TODO: change this to read hdf5
            #profile_matrix_file = open(profile_matrix_path, 'rb')
            #profile_matrix_unpickled = pickle.Unpickler(profile_matrix_file)
            #self.profile_matrix = profile_matrix_unpickled.load()
            pass
        if oma:
            #open oma db object
            ## TODO: unfinished
            #open up taxa Index
            #self.taxtree = ete3.Tree.phylotree('./mastertree.nwk')
            #self.taxaIndex = { n.name:i for i,n in enumerate(self.taxtree.traverse()) }
            #open master tree and taxa index
            with open( './mastertree.pkl', 'rb') as treein:
                self.tree = pickle.loads(treein.read())
            self.tree_string = self.tree.write(format = 1)
            with open( config_utils.datadir + 'taxaIndex.pkl', 'rb') as taxain:
                self.taxaIndex = pickle.loads(taxain.read())
            h5_oma = open_file(config_utils.omadir + 'OmaServer.h5', mode="r")
            self.db_obj = db.Database(h5_oma)
            #open up master tree
            self.treeweights = hashutils.generate_treeweights(self.tree , self.taxaIndex , None, None , None, None)
            self.HAM_PIPELINE = functools.partial(pyhamutils.get_ham_treemap_from_row, tree=self.tree_string )
            self.HASH_PIPELINE = functools.partial(hashutils.row2hash , taxaIndex=self.taxaIndex  , treeweights=self.treeweights , wmg=None )
            self.READ_ORTHO = functools.partial(pyhamutils.get_orthoxml, db_obj=self.db_obj)

    def return_profile_OTF(self, fam, lock = None):
        if type(fam) is str:
            fam = hashutils.hogid2fam(fam)

        if lock is not None:
            lock.acquire()
        ortho_fam = self.READ_ORTHO(fam)
        if lock is not None:
            lock.release()

        tp = self.HAM_PIPELINE([fam, ortho_fam])

        losses = [ self.taxaIndex[n.name]  for n in tp.traverse() if n.lost and n.name in self.taxaIndex  ]
        dupl = [ self.taxaIndex[n.name]  for n in tp.traverse() if n.dupl  and n.name in self.taxaIndex  ]
        presence = [ self.taxaIndex[n.name]  for n in tp.traverse() if n.nbr_genes > 0  and n.name in self.taxaIndex  ]

        indices = dict(zip (['presence', 'loss', 'dup'],[presence,losses,dupl] ) )
        hog_matrix_raw = np.zeros((1, 3*len(self.taxaIndex)))

        for i,event in enumerate(indices):
            if len(indices[event])>0:
                taxindex = np.asarray(indices[event])
                hogindex = np.asarray(indices[event])+i*len(self.taxaIndex)
                hog_matrix_raw[:,hogindex] = 1

        return {fam:{ 'mat':hog_matrix_raw, 'tree':tp} }

    def worker( self,i, inq, retq ):
        print('worker start'+str(i))
        while True:
            input = inq.get()
            if input is None:
                break
            else:
                try:
                    fam,ortho_fam = input
                    tp = self.HAM_PIPELINE([fam, ortho_fam])
                    losses = [ self.taxaIndex[n.name]  for n in tp.traverse() if n.lost and n.name in self.taxaIndex  ]
                    dupl = [ self.taxaIndex[n.name]  for n in tp.traverse() if n.dupl  and n.name in self.taxaIndex  ]
                    presence = [ self.taxaIndex[n.name]  for n in tp.traverse() if n.nbr_genes > 0  and n.name in self.taxaIndex  ]
                    indices = dict(zip (['presence', 'loss', 'dup'],[presence,losses,dupl] ) )
                    hog_matrix_raw = np.zeros((1, 3*len(self.taxaIndex)))
                    for i,event in enumerate(indices):
                        if len(indices[event])>0:
                            taxindex = np.asarray(indices[event])
                            hogindex = np.asarray(indices[event])+i*len(self.taxaIndex)
                            hog_matrix_raw[:,hogindex] = 1
                    retq.put({fam:{ 'mat':hog_matrix_raw, 'tree':tp} })
                except:
                    retq.put({fam:{ 'mat':None, 'tree':None} })


    def retmat_mp(self, traindf , nworkers = 25, chunksize=50  ):
        #fams = [ hashutils.hogid2fam(fam) for fam in fams ]
        def calculate_x(row):
            mat_x1 = row.mat_x
            mat_x2 = row.mat_y

            ret1 = np.zeros(mat_x1.shape)
            ret2 = np.zeros(mat_x2.shape)
            #diff = mat_x1 - mat_x2
            matsum = mat_x1 + mat_x2
            #ret1[np.where(diff != 0 ) ] = -1
            ret2[ np.where(matsum == 2 ) ] = 1
            return list(ret2)


        retq= mp.Queue(-1)
        inq= mp.Queue(-1)
        processes = {}
        mp.log_to_stderr()
        logger = mp.get_logger()
        logger.setLevel(logging.INFO)

        for i in range(nworkers):
            processes[i] = {'time':time.time() , 'process': mp.Process( target = self.worker , args = (i,inq, retq )  ) }
            #processes[i]['process'].daemon = True
            processes[i]['process'].start()

        for batch in range(0, len(traindf) , chunksize ):
            print(batch)

            slicedf = traindf.iloc[batch:batch+chunksize, :]
            fams = list(set(list(slicedf.HogFamA.unique()) + list(slicedf.HogFamB.unique() ) ) )
            total= {}
            for fam in fams:
                orthxml = self.READ_ORTHO(fam)
                if orthxml is not None:
                    inq.put((fam,orthxml))

            done = []
            count = 0
            while len(fams)-1 > count:
                try:
                    data =retq.get(False)
                    count+=1
                    total.update(data)
                except :
                    pass
                time.sleep(.01)

            gc.collect()
            retdf= pd.DataFrame.from_dict( total , orient= 'index')
            slicedf = slicedf.merge( retdf , left_on = 'HogFamA' , right_index = True , how= 'left')
            slicedf = slicedf.merge( retdf , left_on = 'HogFamB' , right_index = True , how= 'left')
            slicedf = slicedf.dropna(subset=['mat_y', 'mat_x'] , how = 'any')

            slicedf['xtrain'] = slicedf.apply( calculate_x , axis = 1)

            X_train = np.vstack( slicedf['xtrain'])
            y_train = slicedf.truth
            print(slicedf)

            yield (X_train, y_train)
        for i in processes:
            inq.put(None)
        for i in processes:
            processes[i]['process'].terminate()



    def hog_query(self, hog_id=None, fam_id=None , k = 100  ):
        """
        Given a hog_id or a fam_id as a query, returns a dictionary containing the results of the LSH.
        :param hog_id: query hog id
        :param fam_id: query fam id
        :return: list containing the results of the LSH for the given query
        """
        if hog_id is not None:
            fam_id = hashutils.hogid2fam(hog_id)
        query_hash = hashutils.fam2hash_hdf5(fam_id, self.hashes_h5 , nsamples=  self.nsamples )

        print(query_hash.hashvalues)

        results = self.lshobj.query(query_hash, k)
        return results


    def hog_query_sorted(self, hog_id=None, fam_id=None , k = 100  ):
        """
        Given a hog_id or a fam_id as a query, returns a dictionary containing the results of the LSH.
        :param hog_id: query hog id
        :param fam_id: query fam id
        :return: list containing the results of the LSH for the given query
        """
        if hog_id is not None:
            fam_id = hashutils.hogid2fam(hog_id)
        query_hash = hashutils.fam2hash_hdf5(fam_id, self.hashes_h5 , nsamples=  self.nsamples )
        results = self.lshobj.query(query_hash, k)
        hogdict = self.pull_hashes(results)

        hogdict = { hog: hogdict[hog].jaccard(query_hash) for hog in hogdict  }


        return hogdict

    def hog_query_OMA(self,hog_id=None, fam_id=None , k = 100 ):
        #construct a profile on the fly
        #rand seed values need to be identical between the construction of the lsh DB and the use of this function
        """
        Untested, Given a hog_id or a fam_id as a query, returns a dictionary containing the results of the LSH.
        Generates the tree profile and hashes on the fly
        :param hog_id: query hog id
        :param fam_id: query fam id
        :return: list containing the results of the LSH for the given query
        """
        if hog_id is not None:
            fam_id = hashutils.hogid2fam(hog_id)
        ortho = self.lshobj.READ_ORTHO(fam)
        tp = self.lshobj.HAM_PIPELINE((fam, ortho))
        hash = self.lshobj.HASH_PIPELINE((fam, tp))
        results = self.lshobj.query(query_hash, k)
        return results

    def pull_hashes(self , hoglist):
        """
        Given a list of hog_ids , returns a dictionary containing their hashes.
        This uses the hdf5 file to get the hashvalues
        :param hog_id: query hog id
        :param fam_id: query fam id
        :return: a dict containing the hash values of the hogs in hoglist
        """
        return { hog: hashutils.fam2hash_hdf5( hashutils.hogid2fam(hog), self.hashes_h5 , nsamples=  self.nsamples) for hog in hoglist}

    def pull_matrows(self,fams):
        """
        given a list of fams return the submatrix containing their profiles

        :return:fams sorted, sparse mat
        """
        return self.profile_matrix[np.asarray(fams),:]


    @staticmethod
    def sort_hashes(query_hash,hashes):
        """
        Given a dict of hogs:hashes, returns a sorted array of hogs and jaccard distances relative to query hog.
        :param query hash: weighted minhash of the query
        :param hashes: a dict of hogs:hashes
        :return: sortedhogs, jaccard
        """
        #sort the hashes by their jaccard relative to query hash
        jaccard=[ query_hash.jaccard(hashes[hog]) for hog in hashes]
        index = np.argsort(jaccard)
        sortedhogs = np.asarry(list(hashes.keys()))[index]
        jaccard= jaccard[index]
        return sortedhogs, jaccard


    @staticmethod
    def allvall_hashes(hashes):
        """
        Given a dict of hogs:hashes, returns generate an all v all jaccard distance matrix.
        :param hashes: a dict of hogs:hashes
        :return: hashmat
        """
        #generate an all v all jaccard distance matrix
        hashmat = np.zeros((len(hashes),len(hashes)))
        for i , hog1 in enumerate(hashes):
            for j, hog2 in enumerate(hashes):
                if i < j :
                    hashmat[i,j]= hashes[hog1].jaccard(hashes[hog2])
        hashmat = hashmat+hashmat.T
        np.fill_diagonal(hashmat, 1)
        return hashmat

    def hog_v_hog(self, hogs):
        """
        give two hogs returns jaccard distance.
        :param hog1 , hog2: str hog id
        :return: jaccard score
        """
        hog1,hog2 = hogs
        #generate an all v all jaccard distance matrix
        hashes = self.pull_hashes([hog1,hog2])
        hashes = list(hashes.values())
        return hashes[0].jaccard(hashes[1])

    def allvall_nx(G,hashes,thresh =None):
        """
        Given a dict of hogs:hashes, returns generate an all v all jaccard distance matrix.
        :param hashes: a dict of hogs:hashes
        :return: hashmat
        """
        #generate an all v all jaccard distance matrix
        hashmat = np.zeros((len(hashes),len(hashes)))

        for i , hog1 in enumerate(hashes):
            for j, hog2 in enumerate(hashes):
                hashmat[i,j]= hashes[hog1].jaccard(hashes[hog2])
        return hashmat

    def iternetwork(seedHOG):
        pass

    def rank_hashes(query_hash,hashes):
        jaccard = []
        sorted = []
        scores = {}
        hogsRanked = np.asarray(list(hashes.keys()))
        for i, hog in enumerate(hashes):
            score = query_hash.jaccard(hashes[hog])
            jaccard.append( score)
            scores[hog] = score
        hogsRanked = list( hogsRanked[ np.argsort(jaccard) ] )
        jaccard = np.sort(jaccard)
        return hogsRanked, jaccard


    def get_vpairs(fam):

        """
        get pairwise distance matrix of OMA all v all
        #not finished
        :param fam: an oma fam
        :return sparesemat: a mat with all taxa in Oma with nonzero entries where this protein is found
        :return densemat: a mat with the taxa covered by the fam
        """
        taxa = self.db_obj.hog_levels_of_fam(fam)
        subtaxindex = { taxon:i for i,taxon in enumerate(taxa)}
        prots = self.db_obj.hog_members_from_hog_id(fam,  'LUCA')
        for prot in prots:
            taxon = prot.ncbi_taxon_id()
            pairs = self.db_obj.get_vpairs(prot)
            for EntryNr1, EntryNr2, RelType , score , distance in list(pairs):
                pass
        return sparsemat , densemat


    def get_submatrix_form_results(self, results):
        res_mat_list = []
        for query, result in results.items():
            res_mat = csr_matrix((len(result), self.profile_matrix.shape[1]))
            for i, r in enumerate(result):
                res_mat[i, :] = self.profile_matrix[r, :]
            res_mat_list.append(res_mat)
        final = np.vstack(res_mat_list)
        return res_mat_list
