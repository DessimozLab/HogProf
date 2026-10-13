import argparse
import os
import pickle
import random
import sys
from collections.abc import Iterable
from pathlib import Path
from typing import overload

import h5py
import numpy as np
import pandas as pd
import tables
from pyoma.browser import db

from hogprof.utils import hashutils, phylo

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

        self.oma = oma
        self.fam2orthoxmlpath = os.path.join(os.path.dirname(lshforestpath), 'fam2orthoxml.csv')
        self._load_hog_fam_mappings()

        ### special operations for sliced subhogs
        if self.slicesubhogs:
            print('sliced subhogs mode')

            # Reconstruct the subhog_id for each row and create a subhog_dict
            subhog_dict = self.id2famsubhog_df.apply(lambda x: f"{x['fam']}_{x['subhog_id']}", axis=1).to_dict()
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
            subhog_to_hashids_dict = {s:h for h, s in self.id2famsubhog_df['subhog_id'].to_dict().items()}
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

        
        self.taxaIndex, self.ReverseTaxaIndex = phylo.generate_taxa_index(self.tree)

        print('DONE with profiler init\n')

    def _load_hog_fam_mappings(self):
        if self.oma and not self.slicesubhogs:
            return
        if os.path.exists(self.fam2orthoxmlpath):
            if not self.slicesubhogs:
                # The CSV index is the builder's HDF5 row ID, not the HOG
                # number or the position of the CSV record.
                mapping = pd.read_csv(self.fam2orthoxmlpath, index_col=0)
                self.hogname_to_fam_dict = {}
                for family, source in mapping['ortho'].items():
                    hog_name = Path(source).stem
                    if hog_name in self.hogname_to_fam_dict:
                        raise ValueError(f"Duplicate rootHOG filename stem: {hog_name}")
                    self.hogname_to_fam_dict[hog_name] = int(family)
                return

            # get a dictionary of subhog ids
            self.id2famsubhog_df = pd.read_csv(self.fam2orthoxmlpath, index_col=None)

            #######################################################################################
            # WARNING:
            # Grouping by `fam` and create a dictionary of indices may not work great for fastOMA where
            # `fam` and `subhog_ids` are completely different
            # Attempt to fix it by Athina:
            # Extract the HOG number (integer) from subhog_id, e.g. "61_HOG:0099866_474" → 99866

            self.id2famsubhog_df["hog_num"] = (
                self.id2famsubhog_df["subhog_id"].str.extract(r"HOG:(\d+)").astype(int)
            )
            fam_dict = (
                self.id2famsubhog_df.groupby("hog_num")
                .apply(lambda x: x.index.tolist())
                .to_dict()
            )
            self.fam_dict = fam_dict

            # Create the reverse dictionary too
            subhogid_to_fam_dict = {
                subhogid: fam
                for fam, subhogid_list in fam_dict.items()
                for subhogid in subhogid_list
            }
            self.subhogid_to_fam_dict = subhogid_to_fam_dict

        elif self.slicesubhogs:
            raise FileNotFoundError(self.fam2orthoxmlpath)
        else:
            # Numeric root IDs do not need a source CSV.
            self.hogname_to_fam_dict = {}
            
    def get_subhogs_of_fam(self, fam_id: int) -> Iterable[int]:
        """
        Returns what sub HOG ID (row numbers in the fam2orthoxml.csv)
        belong to the family of given ID (`fam` column == fam_id)
        """
        if not isinstance(fam_id, int):
            raise TypeError("fam_id must be an integer")

        if not self.slicesubhogs:
            raise ValueError("SubHOGs only supported with --slicesubhogs")

        fam_dict = (
            self.id2famsubhog_df.groupby("fam")
            .apply(lambda x: x.index.tolist())
            .to_dict()
        )
        return fam_dict[fam_id]

    def hogid2fam(self, hog_entry, reverse = False):
        if type(hog_entry) not in [int, np.integer, str]:
            raise TypeError("hog_entry must be integer or string")


        if isinstance(hog_entry, (int, np.integer)) or hog_entry.isdecimal():
            family = int(hog_entry)

        # OMA mode
        elif self.oma:
            if hog_entry.startswith('HOG:'):
                family = int(self.db_obj.parse_hog_id(hog_entry))
            else:
                entry = self.db_obj.entry_by_entry_nr(self.db_obj.id_resolver.resolve(hog_entry))
                family = int(self.db_obj.hog_family(entry=entry))

        # orthoxml mode -- rootHOGs
        elif not self.slicesubhogs and hog_entry in self.hogname_to_fam_dict:
            family = self.hogname_to_fam_dict[hog_entry]

        # orthoxml mode -- subHOGs
        elif self.slicesubhogs and hog_entry in self.subhogname_to_id_dict:
            row = self.subhogname_to_id_dict[hog_entry]
            return [row], {row: hog_entry}

        else:
            raise KeyError(f"Unknown HOG identifier in this OrthoXML index: {hog_entry!r}")

        # return family if the rootHOG mode
        if not self.slicesubhogs:
            return family

        # otherwise return a dict
        rows = list(self.fam_dict[family])
        return rows, {row: self.subhog_dict[row] for row in rows}

    def hog_query(self, *,
                  hog_id: int | None = None,
                  fam_id: int | None = None,
                  k=10):
        """
        Given a hog_id or a fam_id as a query, returns a dictionary containing the results of the LSH.
        :param hog_id: query hog id
        :param fam_id: query fam id
        :return: list containing the results of the LSH for the given query
        """
        if (hog_id is None) == (fam_id is None):
            raise ValueError("Provide either hog_id or fam_id")

        if self.slicesubhogs:
            # resolve fam_id -> list of rows for every subHOG
            subhog_ids = self.get_subhogs_of_fam(fam_id)

            # get subHOG hashes
            query_hashes_dict = {
                i: hashutils.fam2hash_hdf5(i, self.hashes_h5, nsamples=self.nsamples)
                for i in subhog_ids
            }

            # query every subHOG
            return {
                fam: self.lshobj.query(hash, k) for fam, hash in query_hashes_dict.items()
            }
        else:
            fam_id = self.hogid2fam(hog_id if hog_id is not None else fam_id)
            query_hash = hashutils.fam2hash_hdf5(fam_id, self.hashes_h5, nsamples=self.nsamples)
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

        if self.slicesubhogs:
            query_id = hog_id if hog_id is not None else fam_id
            if isinstance(query_id, (int, np.integer)) and query_id not in self.fam_dict:
                return {}, pd.DataFrame()
            family_subhogids_list, _ = self.hogid2fam(query_id)
            fam_id = self.subhogid_to_fam_dict[family_subhogids_list[0]]
        else:
            fam_id = self.hogid2fam(hog_id if hog_id is not None else fam_id)


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

    ##########################################
    # Mini birds
    # fam_id = 133
    # x = p.hog_query(fam_id=fam_id, k=10)
    # print(x)

    # hog_dict, sorted_hogs = p.hog_query_sorted(hog_id=fam_id, k=15000)
    # print(hog_dict)
    # print(sorted_hogs)
    # print("Done")

    # hog_dict, sorted_hogs = p.hog_query_sorted(hog_id="HOG:0005499_1", k=5000)
    # print(hog_dict)

    # print the first and the last HOG from the list of the closest ones
    # for idx in [0, -5, -4, -3, -2, -1]:
    #     hog_id = sorted_hogs[idx]
    #     print(hog_id, hog_dict[hog_id])

    ##########################################
    # Mini birds levels
    fam_id = 102
    x = p.hog_query(fam_id=fam_id, k=10)
    print(x)



if __name__ == '__main__':
    main()
