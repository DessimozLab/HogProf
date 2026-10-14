"""
This test module calls `lshbuilder` on small 
test inputs (orthoXML only). It also does some 
smoke checks if the outputs are as expected.
"""

import csv
import os
import pickle
import shutil
import subprocess
import unittest
import xml.etree.ElementTree as ET
from collections import Counter
from pathlib import Path
from tempfile import TemporaryDirectory

import h5py
import numpy as np
from datasketch import WeightedMinHash
from ete3 import Tree

DATA = Path(__file__).resolve().parent / "data" / "micro_birds"
NS = {"o": "http://orthoXML.org/2011/"}
NPERM = 32
OUTPUT_FILES = {
    "master_tree.corrected.nwk",
    "taxaIndex.pkl",
    "wmg.pkl",
    "hashes.h5",
    "newlshforest.pkl",
    "fam2orthoxml.csv",
}


def family_id(path):
    groups = ET.parse(path).getroot().findall("o:groups/o:orthologGroup", NS)
    if len(groups) != 1:
        raise AssertionError(f"Expected one root family in {path}")
    return groups[0].attrib["id"]


class LshbuilderIntegrationTest(unittest.TestCase):
    def test_root_families(self):
        self.build_and_check(levels=False)

    def test_subhog_families(self):
        self.build_and_check(levels=True)

    def build_and_check(self, *, levels):

        # check that command is present
        executable = shutil.which("lshbuilder")
        self.assertIsNotNone(executable, "Install hogprof to provide lshbuilder")

        # check the input data are there
        inputs = sorted((DATA / "splits").glob("*.orthoxml"))
        self.assertEqual(4, len(inputs), "The complete micro_birds fixture is required")
        
        # check every file has 1 rootHOG
        expected_families = {family_id(path) for path in inputs}
        self.assertEqual(len(inputs), len(expected_families))

        with TemporaryDirectory(prefix="hogprof-ci-") as temp_dir:
            work = Path(temp_dir)
            output = work / "output"
            # Exercise canonical names in root mode and compatibility aliases
            # in levels mode, using the same output checks for both.
            source_flag = "--OrthoGlob" if levels else "--orthoxml-glob"
            tree_flag = "--mastertree" if levels else "--species-tree"
            output_flag = "--outpath" if levels else "--output-dir"
            species_flag = "--specieslim" if levels else "--min-species"
            permutations_flag = "--nperm" if levels else "--num-permutations"
            jobs_flag = "--njobs" if levels else "--jobs"
            command = [
                executable,
                source_flag, str(DATA / "splits" / "*.orthoxml"),
                tree_flag, str(DATA / "species_tree.nwk"),
                output_flag, str(output),
                species_flag, "1",
                permutations_flag, str(NPERM),
                jobs_flag, "2",
            ]
            if levels:
                # Keep even event-free subHOGs so every input family is represented
                command += ["--slicesubhogs", "--eventslim", "-1"]
            result = subprocess.run(
                command, cwd=work, capture_output=True, text=True, timeout=120,
                check=False,
                env=dict(os.environ, PYTHONDONTWRITEBYTECODE="1"),
            )
            # check success
            self.assertEqual(0, result.returncode, result.stdout + result.stderr)

            for name in OUTPUT_FILES:
                with self.subTest(output=name):
                    path = output / name
                    self.assertTrue(path.is_file(), f"Missing {path}\n{result.stderr}")
                    self.assertGreater(path.stat().st_size, 0, f"Empty {path}")

            ############ Check outputs ############
            with (output / "fam2orthoxml.csv").open(newline="") as source:
                rows = list(csv.DictReader(source))

            self.assertTrue(rows, "The family mapping must not be empty")
            mapped_paths = [Path(row["ortho"]).resolve() for row in rows]

            self.assertEqual(set(inputs), set(mapped_paths))
            self.assertEqual(expected_families, {family_id(path) for path in mapped_paths})

            # for the default mode, the number of files = len(mapping)
            if not levels:
                self.assertEqual(Counter(inputs), Counter(mapped_paths))

            # Family numbers are assigned by the CLI; verify each number always
            # points to one input family, including repeated rows in levels mode.
            family_column = "fam" if levels else ""
            families = {}
            for row, path in zip(rows, mapped_paths):
                fam = int(row[family_column])
                self.assertEqual(path, families.setdefault(fam, path))

            self.assertEqual(set(range(len(inputs))), set(families))

            with (output / "taxaIndex.pkl").open("rb") as source:
                taxa = pickle.load(source)

            tree = Tree(
                str(output / "master_tree.corrected.nwk"),
                format=1, quoted_node_names=True,
            )
            self.assertEqual({node.name for node in tree.traverse()}, set(taxa))
            self.assertEqual(len(taxa), len(set(taxa.values())))

            with (output / "wmg.pkl").open("rb") as source:
                generator = pickle.load(source)

            self.assertEqual(NPERM, generator.sample_size)
            self.assertEqual(3 * (max(taxa.values()) + 1), generator.dim)

            with (output / "newlshforest.pkl").open("rb") as source:
                forest = pickle.load(source)

            keys = [
                f"{row['fam']}_{row['subhog_id']}" if levels else row[""]
                for row in rows
            ]
            self.assertEqual(len(keys), len(set(keys)))
            self.assertEqual(set(keys), set(forest.keys))
            self.assertFalse(forest.is_empty())

            with h5py.File(output / "hashes.h5", "r") as hashes:
                self.assertEqual(["NoFilterNoMask"], list(hashes))
                dataset = hashes["NoFilterNoMask"]
                self.assertEqual(2 * NPERM, dataset.shape[1])
                self.assertGreaterEqual(dataset.shape[0], len(keys))

                for index, (row, key) in enumerate(zip(rows, keys)):
                    hash_row = index if levels else int(row[""])

                    # Match the int64 representation used by datasketch and
                    # HogProf's hash reader (HDF5 stores compact int32 values).
                    values = dataset[hash_row].reshape(NPERM, 2).astype(np.int64)
                    self.assertTrue(np.any(values), f"Missing hash for {key}")

                    hashed = WeightedMinHash(seed=generator.seed, hashvalues=values)
                    self.assertIn(key, forest.query(hashed, k=len(keys)))


if __name__ == "__main__":
    unittest.main()
