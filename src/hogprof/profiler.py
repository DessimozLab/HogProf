import argparse
import functools
import os
import pickle
import random
import sys

import h5py
import numpy as np
import pandas as pd
import tables
from pyoma.browser import db

from hogprof.utils import hashutils, phylo, pyhamutils

np.random.seed(0)
random.seed(0)


class Profiler:

    """
    A profiler object allows the user to query the LSH with HOGs and get a list of result HOGs back

    """
    def __init__(self,
                 lshforestpath=None,
                 hashes_h5=None,
                 mat_path=None,
                 oma=False,
                 nsamples=256,
                 mastertree=None,
                 reformat_names=False,
                 swap2taxcode=False,
                 taxfilter=None,
                 taxmask=None,
                 slicesubhogs=False):
        """
        The Profiler class initializes a profiler object for querying the LSH with HOGs and returning a list of result HOGs.

        Attributes:
        lshobj (object): LSH object for querying.
        hashes_h5 (h5py.File): H5 file containing HOGs.
        nsamples (int): Number of samples to use.
        tree (ete3.Tree): Master tree used for generating taxa index.
        tree_string (str): String representation of the master tree.
        taxaIndex (dict): Dictionary mapping taxa names to their indices in the master tree.
        ReverseTaxaIndex (dict): Dictionary mapping indices in the master tree to their corresponding taxa names.
        db_obj (db.Database): OMA database object.
        treeweights (dict): Dictionary containing the tree weight for each taxon.
        READ_ORTHO (callable): Function for reading orthoxml files from OMA.
        HAM_PIPELINE (callable): Function for generating the Annotated tree from a row.
        HASH_PIPELINE (callable): Function for generating the hash from a row.
        slicesubhogs (bool): Whether to slice subhogs.

        Parameters:
        lshforestpath (str, optional): Path to the pickled LSH forest object.
        hashes_h5 (str, optional): Path to the H5 file containing HOGs.
        mat_path (str, optional): Path to the matrix file containing HOGs.
        oma (str, optional): Path to the OMA database.
        tar (str, optional): Path to the tar archive.
        nsamples (int, optional): Number of samples to use. Defaults to 256.
        mastertree (str, optional): Path to the master tree file.
        """
        print('loading lsh')
        with open(lshforestpath, 'rb') as lshpickle:
            self.lshobj = pickle.loads(lshpickle.read())
            print('indexing lsh')
            self.lshobj.index()
        self.hashes_h5 = h5py.File(hashes_h5, mode='r')
        self.slicesubhogs = slicesubhogs
        hamfunction = pyhamutils.get_ham_treemap_from_row
        hashfunction = hashutils.row2hash
        ### special operations for sliced subhogs
        if self.slicesubhogs is True:
            print('sliced subhogs mode')
            #try:
            self.fam2orthoxmlpath = os.path.join(os.path.dirname(lshforestpath), 'fam2orthoxml.csv')
            ### check that it exists
            if not os.path.exists(self.fam2orthoxmlpath):
                print("Could not find fam2orthoxml.csv file. Was HogProf run on OMA?")
                raise FileNotFoundError
            ### adjust functions
            hamfunction = pyhamutils.get_subhog_ham_treemaps_from_row
            hashfunction = hashutils.hash_trees_subhogs
            ### get a dictionary of subhog ids
            id2famsubhog_df = pd.read_csv(self.fam2orthoxmlpath, index_col=None) ### Athina note: here used to be index_col=0 !!!!!!!!!!!!!!!!!!!!!!!!
            # print for fam=15
            print(id2famsubhog_df.head())
            #print(id2famsubhog_df[id2famsubhog_df['fam'] == 15]) 
            # Group by 'fam' and create a dictionary of indices - this may not work great for fastOMA where fam and subhog_ids are completely different!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            #fam_dict = id2famsubhog_df.groupby('fam').apply(lambda x: x.index.tolist()).to_dict()
            #self.fam_dict = fam_dict
            ### attempt to fix it by Athina. # Extract the HOG number (integer) from subhog_id, e.g. "61_HOG:0099866_474" → 99866
            id2famsubhog_df["hog_num"] = id2famsubhog_df["subhog_id"].str.extract(r"HOG:(\d+)").astype(int)
            fam_dict = id2famsubhog_df.groupby("hog_num").apply(lambda x: x.index.tolist()).to_dict()
            self.fam_dict = fam_dict
            # Create the reverse dictionary too
            subhogid_to_fam_dict = {subhogid: fam for fam, subhogid_list in fam_dict.items() for subhogid in subhogid_list}
            self.subhogid_to_fam_dict = subhogid_to_fam_dict
            # Reconstruct the subhog_id for each row and create a subhog_dict
            subhog_dict = id2famsubhog_df.apply(lambda x: f"{x['fam']}_{x['subhog_id']}", axis=1).to_dict()
            self.subhog_dict = subhog_dict
            # Create the reverse dictionary too
            subhogname_to_id_dict = {subhogname: subhogid for subhogid, subhogname in subhog_dict.items()}
            self.subhogname_to_id_dict = subhogname_to_id_dict
            # Connect subhognames to fams
            subhogname_to_fam_dict = {subhog_id: subhog_id.split('_')[0] for subhog_id in subhogname_to_id_dict.keys()}
            self.subhogname_to_fam_dict = subhogname_to_fam_dict
            # Connect subhogids to levels
            subhog_to_level_dict = {subhog_id: subhog_id.split('_')[1] for subhog_id in subhogname_to_id_dict.keys()}
            self.subhog_to_level_dict = subhog_to_level_dict
            # Connect subhogids to hashids (use default row index as hashids)
            subhog_to_hashids_dict = {s:h for h, s in id2famsubhog_df['subhog_id'].to_dict().items()}
            self.subhog_to_hashids_dict = subhog_to_hashids_dict
            #except:
            #   print("Could not find fam2orthoxml.csv file. Was HogProf run on OMA?")
            
            
        print('h5' , self.hashes_h5 , self.hashes_h5.keys())
        self.nsamples = nsamples

        self.tree = phylo.from_file(mastertree)
        print('loaded tree:', self.tree)
        self.tree_string = phylo.to_string(self.tree)
        
        if oma:
            self.reformat_names = reformat_names
            
            self.saving_path = mat_path

            if not os.path.exists(self.saving_path):
                os.makedirs(self.saving_path)
            
            if self.reformat_names:
                raise NotImplementedError("--reformat_names is temporarily disabled")

                # Disabled currently to sort out tree problems one by one
                # in version 0.0.13. Will be either removed as an option or
                # reimplemented later. The original logic is commented below
                #
                # self.tree, self.idmapper = pyhamutils.tree2numerical(self.tree)
                # self.tree_string = self.tree.write(format=1)
                # with open( self.saving_path + 'reformatted_tree.nwk', 'w') as treeout:
                #   treeout.write(self.tree.write(format=0 ))
                # with open( self.saving_path + 'idmapper.pkl', 'wb') as idout:
                #   idout.write( pickle.dumps(self.idmapper))
                # print('reformatted tree')
                # print( self.tree )
                # self.tree_string = self.tree.write(format=1)
                # #remap taxfilter and taxmask
                # if taxfilter:
                #   self.tax_filter = [ self.idmapper[tax] for tax in taxfilter ]
                # if taxmask:
                #   self.tax_mask = self.idmapper[taxmask]
            else:
                self.idmapper = None
                self.tax_filter = taxfilter
                self.tax_mask = taxmask
                self.tree_string = self.tree.write(format=1)
            
            self.taxaIndex, self.ReverseTaxaIndex = phylo.generate_taxa_index(self.tree)
            self.treeweights = hashutils.generate_treeweights(self.tree , self.taxaIndex , None, None )
            self.swap2taxcode = swap2taxcode
            self.tax_filter = None
            self.tax_mask = None
            

            h5_oma = tables.open_file(oma, mode="r")
            self.db_obj = db.Database(h5_oma)
            #self.treeweights = hashutils.generate_treeweights(self.tree , self.taxaIndex , None, None )
            

            self.READ_ORTHO = functools.partial(pyhamutils.get_orthoxml_oma, db_obj=self.db_obj)
            
            self.HAM_PIPELINE = functools.partial(hamfunction, tree=self.tree_string ,  swap_ids=self.swap2taxcode , reformat_names = self.reformat_names , 
                                                  orthoXML_as_string = True, orthomapper = self.idmapper )

            self.HASH_PIPELINE = functools.partial(hashutils.row2hash , taxaIndex=self.taxaIndex  , treeweights=self.treeweights , wmg=None )
        
        self.taxaIndex, self.ReverseTaxaIndex = phylo.generate_taxa_index(self.tree)

        print('DONE with profiler init\n')

    def hogid2fam(self, hog_entry, reverse = False):
        #print("hog_entry", hog_entry)
        if type(hog_entry )== int and not self.slicesubhogs:
            return hog_entry
        elif self.slicesubhogs:
            ### Athina comment: possible redundancy here!!!!!!!!!!!!!!!!!!
            ### if there is a column without label, use it as index 
            #id2famsubhog_df = pd.read_csv(self.fam2orthoxmlpath, index_col=0)
            ### check if index is labeled
            #if id2famsubhog_df.index.name is not None:
            #   id2famsubhog_df = pd.read_csv(self.fam2orthoxmlpath)
            #print(id2famsubhog_df)
            ### fam_dict collects the hashids of all families. it can be used if we want hits for all members of a family
            fam_dict = self.fam_dict#id2famsubhog_df.groupby('fam').apply(lambda x: x.index.tolist()).to_dict()
            #print(fam_dict)
            #subhog_dict = id2famsubhog_df.set_index('subhog_id').to_dict(orient='index')
            if isinstance(hog_entry, int):
                indices = fam_dict[hog_entry]
                if reverse:
                    ### return also the indices to subhogids dict
                    #print(id2famsubhog_df)
                    #subhog_dict = id2famsubhog_df.set_index('subhog_id').to_dict(orient='index')
                    ### for each row of the df with index in indices, reconstruct the subhog id from columns fam and subhog_id
                    subhog_dict = self.subhog_dict#id2famsubhog_df.loc[indices, ['fam', 'subhog_id']].apply(lambda x: f"{x['fam']}_{x['subhog_id']}", axis=1).to_dict()
                    #print(subhog_dict)
                    return indices, subhog_dict
            elif isinstance(hog_entry, str):
                #fam_parts = hog_entry.split('_')
                #fam_int = int(fam_parts[0])
                #subhog_str = '_'.join(fam_parts[1:])
                #indices = id2famsubhog_df[(id2famsubhog_df['fam'] == fam_int) & (id2famsubhog_df['subhog_id'] == subhog_str)].index.tolist()
                indices = [self.subhogname_to_id_dict[ hog_entry ]]
                return indices, {indices[0]:hog_entry}
            else:
                return indices
        else:
            hog_entry = self.db_obj.entry_by_entry_nr(self.db_obj.id_resolver.resolve(hog_entry))
            famnr = int(self.db_obj.hog_family( entry=hog_entry ) )
            return famnr

    def hog_query(self, hog_id=None, fam_id=None, k=10):
        """
        Given a hog_id or a fam_id as a query, returns a dictionary containing the results of the LSH.
        :param hog_id: query hog id
        :param fam_id: query fam id
        :return: list containing the results of the LSH for the given query
        """
        if (hog_id is None) == (fam_id is None):
            raise ValueError("Provide either hog_id or fam_id")

        if hog_id is not None and not self.slicesubhogs:
            fam_id = hog_id

        fam2orthoxmlpath = self.fam2orthoxmlpath if self.slicesubhogs else None
        query_hash = hashutils.fam2hash_hdf5(fam_id, self.hashes_h5,
                                             nsamples=self.nsamples,
                                             fam2orthoxmlpath=fam2orthoxmlpath)
        return self.lshobj.query(query_hash, k)


    def hog_query_sorted(self, hog_id=None, fam_id=None, k=10, fam2orthoxmlpath=None):
        """
        Given a hog_id or a fam_id as a query, returns a dictionary containing the results of the LSH.
        :param hog_id: query hog id
        :param fam_id: query fam id
        :return: list containing the results of the LSH for the given query
        """
        if (hog_id is None) == (fam_id is None):
            raise ValueError("Provide either hog_id or fam_id")

        if hog_id is not None:
            if not self.slicesubhogs:
                fam_id = self.hogid2fam(hog_id)
            else:
                ### if a family is query
                if isinstance(hog_id, int):
                    fam_id = hog_id
                    if not hog_id in self.fam_dict:
                        return {}, pd.DataFrame()
                    family_subhogids_list = self.fam_dict[hog_id]
                    family_subhogs_list = [self.subhog_dict[i] for i in family_subhogids_list]
                ### if a subhog is query
                else:
                    fam_id = int(hog_id.split("_")[0])
                    family_subhogs_list = [hog_id]
                    family_subhogids_list = [self.subhogname_to_id_dict[hog_id]]


        ### case where a family is given as query (including many subfamilies)
        if self.slicesubhogs:
            #print("\nHOG query sorted - Case of sliced subhogs")
            #print(hog_id,'first fam_id', fam_id[0]) # family, [index1, index2, ...]
            query_hashes_dict = {i:hashutils.fam2hash_hdf5(i, self.hashes_h5 , nsamples=  self.nsamples) 
                        for i in family_subhogids_list}
            #print('query_hashes',list(query_hashes_dict.keys())[0],query_hashes_dict[list(query_hashes_dict.keys())[0]]) # index1: hash1, index2: hash2, ...
            ### works without filtering
            results_dict = { fam: self.lshobj.query(query_hash, k) for fam, query_hash in query_hashes_dict.items() }
            #print('results',list(results_dict.keys())[0],'first result hogid',results_dict[list(results_dict.keys())[0]][0]) # index1: results1, index2: results2, ...
            ### the line below is not optimised!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            hogdict = { fam: self.pull_hashes(results, fam2orthoxmlpath=fam2orthoxmlpath) for fam, results in results_dict.items() }
            #print('hogdict',list(hogdict.keys())[0],'first result: 0_0_0',hogdict[list(hogdict.keys())[0]]['0_0_0']) # {index1:{hitsubhog1:hash1, hitsubhog2:hash2, ...}, index2:{subhog1:hash1, subhog2:hash2, ...}, ...}             
            hogdict = { fam: { subhog: hogdict[fam][subhog].jaccard(query_hashes_dict[fam]) for subhog in hogdict[fam] } for fam in hogdict }
            #print('hogdict',list(hogdict.keys())[0],'first result: 0_0_0',hogdict[list(hogdict.keys())[0]]['0_0_0']) # {index1: jaccard1, index2: jaccard2, ...}

            ### making df instead of sorted list to include family info
            data = []
            for index, hit_dict in hogdict.items():
                for subhoghit, jaccard in hit_dict.items():
                    data.append({
                        'query_subhogid': self.subhog_dict[index],
                        'hit_subhogid': subhoghit,
                        'jaccard': jaccard
                    })
            subhogsdf = pd.DataFrame(data, columns=['query_subhogid', 'hit_subhogid', 'jaccard'])
            subhogsdf['query_family'] = fam_id
            #print(subhogsdf.head())
            # Add query_family and query_level columns using the dictionaries
            subhogsdf['query_level'] = subhogsdf['query_subhogid'].map(self.subhog_to_level_dict)
            # Add hit_family and hit_level columns using the dictionaries
            subhogsdf['hit_family'] = subhogsdf['hit_subhogid'].map(self.subhogname_to_fam_dict)
            subhogsdf['hit_level'] = subhogsdf['hit_subhogid'].map(self.subhog_to_level_dict)
            subhogsdf = subhogsdf.sort_values(['jaccard', 'query_subhogid', 'hit_family', 'hit_level'], ascending=False)
            return hogdict, subhogsdf
        else:
            query_hash = hashutils.fam2hash_hdf5(fam_id, self.hashes_h5, nsamples=self.nsamples)
            closest_hits = [int(x) for x in self.lshobj.query(query_hash, k)]

            # pull family hashes
            fam2hash = lambda fam_id: hashutils.fam2hash_hdf5(fam_id, self.hashes_h5,
                                                              nsamples=self.nsamples)
            hit2hash = { hit_id: fam2hash(hit_id) for hit_id in closest_hits }

            # compute jaccard similarity to every hit
            hit2dist = {
                hit_id: hit_hash.jaccard(query_hash)
                for hit_id, hit_hash in hit2hash.items()
            }
            sorted_hogs = sorted(hit2dist, key=hit2dist.get, reverse=True)
            return hit2dist, sorted_hogs

    def pull_hashes(self, hoglist, fam2orthoxmlpath = None):

        """
        Given a list of hog_ids , returns a dictionary containing their hashes.
        This uses the hdf5 file to get the hashvalues
        :param hog_id: query hog id
        :param fam_id: query fam id
        :return: a dict containing the hash values of the hogs in hoglist
        """
        #print('pulling hashes')
        ### hogid2fam returns the index of the family in the h5 file for normal cases
        ### for sliced subhogs, it returns the list of indices of the subhogs in the h5 file
        ### returns indices, subhog_dict
        if not self.slicesubhogs:
            return { entry: hashutils.fam2hash_hdf5( self.hogid2fam(entry), self.hashes_h5 , nsamples=  self.nsamples) for entry in hoglist}
        else:
            result_dict = {}
            for entry in hoglist:
                indices, subhog_dict = self.hogid2fam(entry)
                result_dict[entry] = [hashutils.fam2hash_hdf5(idx, self.hashes_h5, nsamples=self.nsamples) for idx in indices][0]
            return result_dict


    @staticmethod
    def sort_hashes(query_hash, hashes):
        """
        Given a dict of hogs:hashes, returns a sorted array of hogs and jaccard distances relative to query hog.
        :param query hash: weighted minhash of the query
        :param hashes: a dict of hogs:hashes
        :return: sortedhogs, jaccard
        """
        # sort the hashes by their jaccard relative to query hash
        jaccard = [query_hash.jaccard(hashes[hog]) for hog in hashes]
        index = np.argsort(jaccard)
        sorted_hogs = np.asarray(list(hashes.keys()))[index]
        jaccard = jaccard[index]
        return sorted_hogs, jaccard


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


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--lshforestpath", help="Path to the newlshforest.plk file", type=str
    )
    parser.add_argument("--hashes_h5", help="Path to the hashes.h5 file", type=str)
    parser.add_argument(
        "--mat_path", help="Path to the fam2orthoxml.csv file", type=str
    )
    parser.add_argument("--OMA", help="Path to OmaServer.h5 file", type=str)
    parser.add_argument("--nsamples", help="no of samples", type=int, default=256)
    parser.add_argument(
        "--slicesubhogs", help="Use slicesubhogs mode", action="store_true"
    )
    parser.add_argument(
        "--mastertree",
        help="master taxonomic tree. nodes should correspond to orthoxml",
        type=str,
    )
    args = vars(parser.parse_args(sys.argv[1:]))
    print("\nReading arguments")

    if args["slicesubhogs"] == "True":
        args["slicesubhogs"] = True
    elif args["slicesubhogs"] == "False":
        args["slicesubhogs"] = False

    p = Profiler(
        lshforestpath=args["lshforestpath"],
        hashes_h5=args["hashes_h5"],
        mat_path=args["mat_path"],
        oma=args["OMA"],
        nsamples=args["nsamples"],
        slicesubhogs=args["slicesubhogs"],
        mastertree=args["mastertree"],
    )
    print("\nProfiler object created")

    #fam_id = 133
    #x = p.hog_query(fam_id=fam_id, k=10)
    #print(x)

    #hog_dict, sorted_hogs = p.hog_query_sorted(hog_id=fam_id, k=15000)
    #print(hog_dict)
    #print(sorted_hogs)
    #print("Done")
    
    #hogdict, sortedhogs = p.hog_query_sorted( hog_id= "XP_015681225" , k = 20 )
    #print(sortedhogs)
    #print()
    #if sortedhogs['hit_subhogid'].str.contains('0_4_0').any(): ## for local was 0_0_0, for curnagl was 0_4_0
    #   print('got hit!\n')
    #else:
    #   print('Warning! Did not find itself!\n')


if __name__ == '__main__':
    main()
