"""HogEntryIndex lookups."""

import unittest
from pathlib import Path
from tempfile import TemporaryDirectory

import numpy as np
import polars as pl

from hogprof.metadata import HogEntry, HogEntryIndex


def write_mapping(path: Path, rows: list[dict]) -> Path:
    pl.DataFrame(rows, schema=HogEntryIndex.SCHEMA).write_csv(path)
    return path


class HogEntryIndexTest(unittest.TestCase):
    def setUp(self):
        directory = TemporaryDirectory()
        self.addCleanup(directory.cleanup)
        self.dir = Path(directory.name)

        # Keys out of order and with a gap, like a rootHOG-mode index; one
        # HOG id (HOG:1.1a) defined at two levels; levels interleaved
        self.rows = [
            {"hog_key": 7, "fam_id": "HOG:1", "hog_id": "HOG:1", "level": "root", "source": "a"},
            {"hog_key": 2, "fam_id": "HOG:1", "hog_id": "HOG:1.1a", "level": "Aves", "source": "a"},
            {"hog_key": 5, "fam_id": "HOG:1", "hog_id": "HOG:1.1a", "level": "Neognathae", "source": "a"},
            {"hog_key": 0, "fam_id": "HOG:2", "hog_id": "HOG:2", "level": "root", "source": "b"},
            {"hog_key": 3, "fam_id": "HOG:2", "hog_id": "HOG:2.1b", "level": "Aves", "source": "b"},
        ]
        self.index = HogEntryIndex(write_mapping(self.dir / "map.csv", self.rows))

    def test_get_hog_by_key(self):
        self.assertFalse(self.index._key_is_row)
        for row in self.rows:
            self.assertEqual(HogEntry(**row), self.index.get_hog_by_key(row["hog_key"]))

    def test_missing_and_invalid_keys(self):
        for key in (1, 4, 6, 8, -1):
            with self.subTest(key=key), self.assertRaises(ValueError):
                self.index.get_hog_by_key(key)
        with self.assertRaises(TypeError):
            self.index.get_hog_by_key(np.int64(2))

    def test_get_hog_by_hog_id(self):
        hogs = self.index.get_hog_by_hog_id("HOG:1.1a")
        self.assertEqual(["Aves", "Neognathae"], [h.level for h in hogs])
        self.assertEqual([], self.index.get_hog_by_hog_id("HOG:404"))

    def test_levels(self):
        self.assertEqual(["Aves", "Neognathae", "root"], self.index.levels())
        self.assertEqual(5, len(self.index))

        # file order is kept
        self.assertEqual(["HOG:1.1a", "HOG:2.1b"], [h.hog_id for h in self.index.search_by_level("Aves")])
        np.testing.assert_array_equal([7, 0], self.index.keys_at_level("root"))
        self.assertEqual([7, 0], self.index.get_level("root")["hog_key"].to_list())

        self.assertEqual([], self.index.search_by_level("Mammalia"))
        self.assertEqual(0, self.index.keys_at_level("Mammalia").size)
        self.assertEqual(0, self.index.get_level("Mammalia").height)

    def test_rejects_duplicate_keys(self):
        rows = [*self.rows, {**self.rows[0], "hog_id": "HOG:dup"}]
        with self.assertRaisesRegex(ValueError, "duplicate hog_key 7"):
            HogEntryIndex(write_mapping(self.dir / "dup.csv", rows))


class KeyIsRowTest(unittest.TestCase):
    """subHOG indexes: hog_key == row number, looked up without a search"""

    def setUp(self):
        directory = TemporaryDirectory()
        self.addCleanup(directory.cleanup)
        self.rows = [
            {"hog_key": i, "fam_id": f"HOG:{i // 2}", "hog_id": f"HOG:{i // 2}.{i}",
             "level": "root" if i % 2 else "Aves", "source": ""}
            for i in range(6)
        ]
        self.index = HogEntryIndex(write_mapping(Path(directory.name) / "map.csv", self.rows))

    def test_key_is_row(self):
        self.assertTrue(self.index._key_is_row)
        self.assertIsNone(self.index._sorted_keys)
        for row in self.rows:
            self.assertEqual(row["hog_id"], self.index.get_hog_by_key(row["hog_key"]).hog_id)
        np.testing.assert_array_equal([0, 2, 4], self.index.keys_at_level("Aves"))

    def test_out_of_range_keys(self):
        for key in (-1, 6, 100):
            with self.subTest(key=key), self.assertRaises(ValueError):
                self.index.get_hog_by_key(key)

    def test_unordered_keys_use_search(self):
        # same keys, not in row order
        directory = TemporaryDirectory()
        self.addCleanup(directory.cleanup)
        index = HogEntryIndex(write_mapping(Path(directory.name) / "map.csv", self.rows[::-1]))
        self.assertFalse(index._key_is_row)
        for row in self.rows:
            self.assertEqual(row["hog_id"], index.get_hog_by_key(row["hog_key"]).hog_id)


if __name__ == "__main__":
    unittest.main()
