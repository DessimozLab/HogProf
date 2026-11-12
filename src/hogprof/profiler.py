import logging
import os
import pickle
import random
from collections.abc import Iterable
from pathlib import Path

import numpy as np

from hogprof import hashstore
from hogprof.artifacts import FOREST, HASHSTORE, MAPPING
from hogprof.hashstore import HashStore
from hogprof.manifest import IndexManifest
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
                 nsamples=256,
                 slicesubhogs=False,
                 hash_seed=hashutils.HASH_SEED,
                 hashes_dataset=None):
        """
        The Profiler class initializes a profiler object for querying the
        LSH with HOGs and returning a list of result HOGs.

        Prefer `Profiler.from_index_dir`, which reads these parameters
        from the index manifest.

        lshforestpath (str, optional): Path to the pickled LSH forest object.
        hashes_h5 (Path): hashes.h5 of the index.
        nsamples (int, optional): Number of permutations of the index. Defaults to 256.
        slicesubhogs (bool): Whether the index has subHOG profiles.
        hash_seed (int): Seed of the hash generator used to build the index.
        hashes_dataset (str, optional): Dataset name in `hashes_h5`; the first one if None.
        """
        logger.debug('Loading LSH')
        with open(lshforestpath, 'rb') as lshpickle:
            self.forest = pickle.loads(lshpickle.read())
            logger.debug('Indexing LSH')

        self.hashes = HashStore(hashes_h5, num_permutations=nsamples,
                                dataset=hashes_dataset, seed=hash_seed)
        self.slice_subhogs = slicesubhogs

        self.fam2orthoxmlpath = os.path.join(os.path.dirname(lshforestpath), MAPPING)

        self.hog_index = HogEntryIndex(Path(self.fam2orthoxmlpath))

        logger.debug("Hashes: %s, dataset %s", self.hashes.path, self.hashes.dataset_name)
        self.nsamples = nsamples

        logger.debug('Profiler init over')

    @classmethod
    def from_index_dir(cls, index_dir) -> "Profiler":
        """Open an index built by `hogprof build`, configured from its manifest"""
        index_dir = Path(index_dir)
        manifest = IndexManifest.read(index_dir)
        return cls(
            lshforestpath=index_dir / FOREST,
            hashes_h5=index_dir / HASHSTORE,
            nsamples=manifest.num_permutations,
            slicesubhogs=manifest.slice_subhogs,
            hash_seed=manifest.hash_seed,
            hashes_dataset=manifest.hashes_dataset,
        )

    def close(self) -> None:
        self.hashes.close()

    def __enter__(self) -> "Profiler":
        return self

    def __exit__(self, *exc) -> None:
        self.close()

    def get_subhogs_of_fam(self, fam_id: int) -> Iterable[int]:
        """
        Returns what sub HOG ID (row numbers in the fam2orthoxml.csv)
        belong to the family of given ID (`fam` column == fam_id)
        """
        if not isinstance(fam_id, str):
            raise TypeError("fam_id must be a string")

        if not self.slice_subhogs:
            raise ValueError("SubHOGs only supported with --slicesubhogs")

        fam_dict = (
            self.id2famsubhog_df.groupby("fam")
            .apply(lambda x: x.index.tolist())
            .to_dict()
        )
        return fam_dict[fam_id]


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

        query_packed = self.hashes.get_one(hog.hog_key)
        query_hash = hashstore.to_minhash(query_packed, seed=self.hashes.seed)

        # FIXME: there must be a better way of doing
        # this, e.g. with splitting the forest into sets
        # of forests, one per level, and querying exactly k
        # in the required level-defined forest.
        # Current solution: for levels, let's query twice more
        # hits, since they might be filtered out further
        # if they are from a wrong level.
        queried_k = 2 * k if level else k

        hits = []
        for hit_key in self.forest.query(query_hash, queried_k):
            hit = self.hog_index.get_hog_by_key(int(hit_key))

            # If --level is provided, we need to make sure that
            # neighbors are defined at the same level as the query.
            if level and level != hit.level:
                continue

            hits.append(hit)
            if len(hits) >= k:
                break

        # one bulk read and comparison for all hits
        jaccards = hashstore.jaccard(query_packed, self.hashes.get([h.hog_key for h in hits]))

        rows = [
            HogQueryHit(
                query_hog_id=hog_name,
                hit_hog_id=hit.hog_id,
                hit_level=hit.level,
                jaccard=float(j),
            )
            for hit, j in zip(hits, jaccards)
        ]
        return HogQueryResult.from_hits(rows)


    def query_hog_significant(self,
                              hog_name: str,
                              level : str | None = None,
                              threshold: float = 0.05,
                              max_candidates: int = 100) -> HogQueryResult:
        """
        Find all significantly close HOGs to the one with `hog_id`.
        """
        raise NotImplementedError()


def run_search(index_dir: Path,
               hog_id: str,
               level: str = None,
               k: int = 10,
               ):
    p = Profiler.from_index_dir(index_dir)
    return p.query_hog(hog_id, level=level, k=k).hits