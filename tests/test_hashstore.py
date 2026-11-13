"""HashStore reads and Jaccard estimates match datasketch."""

import itertools
import unittest
from pathlib import Path
from tempfile import TemporaryDirectory

import h5py
import numpy as np
import polars as pl
from datasketch import WeightedMinHashGenerator

from hogprof import hashstore
from hogprof.artifacts import HASHSTORE, MAPPING
from hogprof.hashstore import HashStore
from hogprof.manifest import IndexManifest
from hogprof.metadata import HogEntryIndex

NPERM = 64
DATASET = "NoFilterNoMask"


class HashStoreTest(unittest.TestCase):
    def setUp(self):
        directory = TemporaryDirectory()
        self.addCleanup(directory.cleanup)
        self.dir = Path(directory.name)

        # Binary profiles over 30 columns with heavy overlap, so that
        # pairs cover a range of Jaccard values, including identical ones
        rng = np.random.default_rng(0)
        self.wmg = WeightedMinHashGenerator(30, sample_size=NPERM, seed=1)
        base = rng.random(30) < 0.5
        profiles = [base ^ (rng.random(30) < p) for p in np.linspace(0, 0.5, 12)]
        profiles = [np.where(p.any(), p, base) for p in profiles]
        self.minhashes = [self.wmg.minhash(p.astype(float)) for p in profiles]

        # Store like the saver: int32, raveled (k, t) pairs, with
        # unused zero rows (row 3 is a gap, plus trailing chunk padding)
        self.keys = [0, 1, 2, *range(4, 4 + len(self.minhashes) - 3)]
        with h5py.File(self.dir / HASHSTORE, "w") as h5:
            ds = h5.create_dataset(DATASET, (max(self.keys) + 10, 2 * NPERM), dtype="int32")
            for key, mh in zip(self.keys, self.minhashes):
                ds[key, :] = mh.hashvalues.ravel()

        self.store = HashStore(self.dir / HASHSTORE, num_permutations=NPERM)
        self.addCleanup(self.store.close)

    def test_jaccard_matches_datasketch(self):
        packed = self.store.get(self.keys)
        for (i, a), (j, b) in itertools.combinations(enumerate(self.minhashes), 2):
            expected = a.jaccard(b)
            self.assertEqual(expected, hashstore.batch_jaccard(packed[i], packed[j:j + 1])[0])

        matrix = hashstore.jaccard_matrix(packed)
        expected = np.array([[a.jaccard(b) for b in self.minhashes] for a in self.minhashes])
        np.testing.assert_array_equal(expected, matrix)

        # the test data must not be trivial
        self.assertGreater(len(np.unique(matrix)), 5)

    def test_jaccard_matrix_chunks(self):
        packed = self.store.get(self.keys)
        np.testing.assert_array_equal(
            hashstore.jaccard_matrix(packed, packed[:5]),
            hashstore.jaccard_matrix(packed, packed[:5], chunk_bytes=1),
        )

    def test_minhash_round_trip(self):
        for key, mh in zip(self.keys, self.minhashes):
            restored = self.store.get_minhash(key)
            np.testing.assert_array_equal(mh.hashvalues, restored.hashvalues)
            self.assertEqual(1.0, restored.jaccard(mh))
            np.testing.assert_array_equal(self.store.get_one(key), hashstore.from_minhash(mh))

    def test_get_any_order_and_duplicates(self):
        keys = [self.keys[5], self.keys[0], self.keys[5], self.keys[2]]
        packed = self.store.get(keys)
        self.assertEqual((4, NPERM), packed.shape)
        for row, key in zip(packed, keys):
            np.testing.assert_array_equal(self.store.get_one(key), row)
        self.assertEqual((0, NPERM), self.store.get([]).shape)

    def test_sparse_and_dense_reads_agree(self):
        # dense keys are read as one slice, sparse ones by fancy indexing
        with HashStore(self.dir / HASHSTORE, num_permutations=NPERM, preload=True) as loaded:
            for keys in ([self.keys[0], self.keys[-1]], self.keys[2:7], self.keys[::-1]):
                with self.subTest(keys=keys):
                    np.testing.assert_array_equal(loaded.get(keys), self.store.get(keys))

    def test_preload_is_equivalent(self):
        with HashStore(self.dir / HASHSTORE, num_permutations=NPERM, preload=True) as loaded:
            self.assertTrue(loaded.preloaded)
            keys = self.keys[::-1]
            np.testing.assert_array_equal(self.store.get(keys), loaded.get(keys))

    def test_invalid_keys(self):
        for keys in ([-1], [self.store.num_rows]):
            with self.subTest(keys=keys), self.assertRaises(KeyError):
                self.store.get(keys)
        with self.assertRaises(TypeError):
            self.store.get([0.5])

    def test_rejects_wrong_layout(self):
        with self.assertRaisesRegex(ValueError, "expected 256"):
            HashStore(self.dir / HASHSTORE, num_permutations=128)
        with self.assertRaisesRegex(ValueError, "no dataset"):
            HashStore(self.dir / HASHSTORE, num_permutations=NPERM, dataset="missing")

    def test_from_index_dir(self):
        IndexManifest(num_permutations=NPERM, hash_seed=1, hashes_dataset=DATASET,
                      slice_subhogs=False, num_hogs=len(self.keys)).write(self.dir)
        with HashStore.from_index_dir(self.dir) as store:
            self.assertEqual(DATASET, store.dataset_name)
            np.testing.assert_array_equal(self.store.get(self.keys), store.get(self.keys))


class ReorderByLevelTest(unittest.TestCase):
    def setUp(self):
        directory = TemporaryDirectory()
        self.addCleanup(directory.cleanup)
        self.dir = Path(directory.name)
        self.hashes = self.dir / HASHSTORE
        self.mapping = self.dir / MAPPING

        # Like the saver: levels interleaved family by family, keys 0..n-1,
        # plus unused trailing rows. Each row's values encode its hog_id.
        self.levels = ["root", "Aves", "Neognathae", "Aves", "root", "Neognathae", "Aves", "root"]
        self.rows = [
            {"hog_key": key, "fam_id": f"HOG:{key // 3}", "hog_id": f"HOG:{key // 3}.{key}",
             "level": level, "source": ""}
            for key, level in enumerate(self.levels)
        ]
        self.values = {row["hog_id"]: np.full(2 * NPERM, 100 + row["hog_key"], dtype=np.int32)
                       for row in self.rows}
        with h5py.File(self.hashes, "w") as h5:
            ds = h5.create_dataset(DATASET, (len(self.rows) + 5, 2 * NPERM),
                                   maxshape=(None, 2 * NPERM), dtype="int32")
            for row in self.rows:
                ds[row["hog_key"]] = self.values[row["hog_id"]]
        pl.DataFrame(self.rows).write_csv(self.mapping)

    def test_levels_are_contiguous(self):
        hashstore.reorder_by_level(self.hashes, self.mapping, DATASET)

        index = HogEntryIndex(self.mapping)
        mapping = pl.read_csv(self.mapping, schema_overrides=HogEntryIndex.SCHEMA)
        self.assertEqual(list(range(len(self.rows))), mapping["hog_key"].to_list())
        self.assertEqual(sorted(self.levels), mapping["level"].to_list())

        # same HOGs, and every HOG still points to its own hash
        self.assertEqual({r["hog_id"] for r in self.rows}, set(mapping["hog_id"]))
        with h5py.File(self.hashes, "r") as h5:
            ds = h5[DATASET]
            self.assertEqual((len(self.rows), 2 * NPERM), ds.shape)
            self.assertIsNone(ds.chunks)
            for key, hog_id in zip(mapping["hog_key"], mapping["hog_id"]):
                np.testing.assert_array_equal(self.values[hog_id], ds[key])

        for level in set(self.levels):
            keys = [h.hog_key for h in index.search_by_level(level)]
            self.assertEqual(list(range(keys[0], keys[0] + len(keys))), keys)
            # within a level, the original (family) order is kept
            original = [r["hog_id"] for r in self.rows if r["level"] == level]
            self.assertEqual(original, [h.hog_id for h in index.search_by_level(level)])

    def test_rejects_keys_beyond_hashes(self):
        rows = [*self.rows, {**self.rows[0], "hog_key": 100, "hog_id": "HOG:x"}]
        pl.DataFrame(rows).write_csv(self.mapping)
        with self.assertRaisesRegex(ValueError, "beyond"):
            hashstore.reorder_by_level(self.hashes, self.mapping, DATASET)


if __name__ == "__main__":
    unittest.main()
