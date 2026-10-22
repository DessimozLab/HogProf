import logging
import os
import pickle
import random
from collections.abc import Iterable
from pathlib import Path

import h5py
import numpy as np
import tables
from pyoma.browser import db

from hogprof.metadata import HogEntry, HogEntryIndex, HogQueryHit, HogQueryResult
from hogprof.utils import hashutils

np.random.seed(0)
random.seed(0)
logger = logging.getLogger(__name__)


class Profiler:
    """
    A profiler object allows the user to query the LSH with HOGs and get a list of result HOGs back
    """
    def __init__(self,
                 lshforestpath=None,
                 hashes_h5=None,
                 oma=False,
                 nsamples=256,
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
        logger.debug('Loading LSH')
        with open(lshforestpath, 'rb') as lshpickle:
            self.forest = pickle.loads(lshpickle.read())
            logger.debug('Indexing LSH')
            self.forest.index()

        self.hashes_h5 = h5py.File(hashes_h5, mode='r')
        self.slicesubhogs = slicesubhogs

        self.oma = oma
        self.fam2orthoxmlpath = f"{os.path.dirname(lshforestpath)}/fam2orthoxml.csv"

        self.hog_index = HogEntryIndex(Path(self.fam2orthoxmlpath))

        logger.debug("H5" + str(self.hashes_h5) + " " + str(self.hashes_h5.keys()))
        self.nsamples = nsamples
        
        if oma:
            h5_oma = tables.open_file(oma, mode="r")
            self.db_obj = db.Database(h5_oma)

        logger.debug('Profiler init over')

    def get_subhogs_of_fam(self, fam_id: int) -> Iterable[int]:
        """
        Returns what sub HOG ID (row numbers in the fam2orthoxml.csv)
        belong to the family of given ID (`fam` column == fam_id)
        """
        if not isinstance(fam_id, str):
            raise TypeError("fam_id must be a string")

        if not self.slicesubhogs:
            raise ValueError("SubHOGs only supported with --slicesubhogs")

        fam_dict = (
            self.id2famsubhog_df.groupby("fam")
            .apply(lambda x: x.index.tolist())
            .to_dict()
        )
        return fam_dict[fam_id]

    def hogid2fam(self, hog_entry):
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
        elif not self.slicesubhogs and hog_entry in self.hog2fam:
            family = self.hog2fam[hog_entry]

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


    def _get_hog_of_level(self, hog_entries: list[HogEntry], level: str) -> HogEntry:
        for h in hog_entries:
            if h.level == level:
                return h

        raise ValueError(f"Cannot find the HOG at level {level}")


    def query_hog(self, hog_name: str, level=None, k: int = 10) -> HogQueryResult:
        """
        Query the LSH index directly. `hog_id` MUST be a forest key
        """
        if not isinstance(hog_name, str) or not hog_name.startswith('HOG:'):
            raise TypeError("`hog_name` must be a string of format 'HOG:{number}'")

        hogs = self.hog_index.get_hog_by_hog_id(hog_name)
        if len(hogs) == 0:
            raise ValueError(f"Unknown HOG identifier: {hog_name!r}")

        if len(hogs) > 1:
            if not level:
                levels = [h.level for h in hogs]
                raise ValueError(f"{hog_name} exists at different levels: {levels}. Please indicate a level")

            hog = self._get_hog_of_level(hogs, level)
        else:
            hog = hogs[0]

        query_hash = hashutils.fam2hash_hdf5(hog.hog_key, self.hashes_h5, nsamples=self.nsamples)

        rows = []
        hit_keys = self.forest.query(query_hash, k)
        for hit_key in hit_keys:
            hit_key = int(hit_key)

            # get hit metadata
            hit = self.hog_index.get_hog_by_key(hit_key)

            # compute jaccard
            hit_hash = hashutils.fam2hash_hdf5(hit_key, self.hashes_h5, nsamples=self.nsamples)
            jaccard = query_hash.jaccard(hit_hash)

            # report
            rows.append(
                HogQueryHit(
                    query_hog_id=hog_name,
                    hit_hog_id=hit.hog_id,
                    hit_level=hit.level,
                    jaccard=jaccard,
                )
            )

        hqr = HogQueryResult.from_hits(rows)
        hqr.hits.sort("jaccard", descending=True)
        return hqr


    def query_fam(self, fam_id: str, k=10):
        """
        Query HOGs of a given family.

        Default mode: only the rootHOG will be queried, and
        `fam_id` MUST be a valid forest key. O

        In the --slicesubhogs mode: all subHOGs of the given
        family are queried.

        Returns a dict of type
        {
            HOG_id : [ list of closest neibhbor IDs ]
        }
        """
        raise NotImplementedError()

        if self.slicesubhogs:
            # resolve fam_id -> list of forest key for every subHOG
            subhog_keys = self.get_subhogs_of_fam(fam_id)

            # get subHOG hashes
            query_hashes_dict = {
                key: hashutils.fam2hash_hdf5(key, self.hashes_h5, nsamples=self.nsamples)
                for key in subhog_keys
            }
        else:
            query_hashes_dict = {fam_id: hashutils.fam2hash_hdf5(fam_id, self.hashes_h5, nsamples=self.nsamples)}

        # query every subHOG
        return {
            forest_key: self.forest.query(hash, k)
            for forest_key, hash in query_hashes_dict.items()
        }

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
        logger.warning("`hog_query` is deprecated, use `query_hog` or `query_fam` instead")

        if hog_id:
            return self.query_hog(hog_id, k=k)
        elif fam_id:
            return self.query_fam(fam_id, k=k)
        else:
            raise ValueError("Provide either hog_id or fam_id")


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


def run_search(index_dir: Path,
               hog_id: str,
               level: str = None,
               k: int = 10,
               oma=None,
               nsamples=256,
               slice_subhogs=False,
               ):

    index_dir = Path(index_dir)
    p = Profiler(
        lshforestpath=index_dir / "newlshforest.pkl",
        hashes_h5=index_dir / "hashes.h5",
        oma=oma,
        nsamples=nsamples,
        slicesubhogs=slice_subhogs,
    )
    return p.query_hog(hog_id, level=level, k=k).hits