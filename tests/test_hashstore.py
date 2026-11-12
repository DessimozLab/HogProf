"""HashStore reads and Jaccard estimates match datasketch."""

import itertools
import unittest
from pathlib import Path
from tempfile import TemporaryDirectory

import h5py
import numpy as np
from datasketch import WeightedMinHashGenerator

from hogprof import hashstore
from hogprof.artifacts import HASHSTORE
from hogprof.hashstore import HashStore
from hogprof.manifest import IndexManifest

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


if __name__ == "__main__":
    unittest.main()
