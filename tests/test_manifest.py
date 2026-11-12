"""Index manifest round trip and validation."""

import json
import unittest
from pathlib import Path
from tempfile import TemporaryDirectory

from hogprof.artifacts import MANIFEST
from hogprof.manifest import FORMAT_VERSION, IndexManifest


class IndexManifestTest(unittest.TestCase):
    def setUp(self):
        directory = TemporaryDirectory()
        self.addCleanup(directory.cleanup)
        self.dir = Path(directory.name)
        self.manifest = IndexManifest(
            num_permutations=32, hash_seed=1, hashes_dataset="NoFilterNoMask",
            slice_subhogs=True, num_hogs=7, build={"min_events": 0},
        )

    def test_round_trip(self):
        self.manifest.write(self.dir)
        read = IndexManifest.read(self.dir)
        self.assertEqual(32, read.num_permutations)
        self.assertTrue(read.slice_subhogs)
        self.assertEqual({"min_events": 0}, read.build)
        self.assertTrue(read.created)
        self.assertTrue(read.hogprof_version)
        self.assertFalse((self.dir / (MANIFEST + ".tmp")).exists())

    def test_missing_manifest(self):
        with self.assertRaisesRegex(FileNotFoundError, "rebuild"):
            IndexManifest.read(self.dir)

    def test_rejects_other_format_versions(self):
        self.manifest.write(self.dir)
        path = self.dir / MANIFEST
        data = json.loads(path.read_text())
        data["format_version"] = FORMAT_VERSION + 1
        path.write_text(json.dumps(data))
        with self.assertRaisesRegex(ValueError, "format version"):
            IndexManifest.read(self.dir)

    def test_rejects_missing_fields_and_ignores_unknown(self):
        self.manifest.write(self.dir)
        path = self.dir / MANIFEST
        data = json.loads(path.read_text())
        data["future_field"] = "ignored"
        path.write_text(json.dumps(data))
        self.assertEqual(7, IndexManifest.read(self.dir).num_hogs)

        del data["num_permutations"]
        path.write_text(json.dumps(data))
        with self.assertRaisesRegex(ValueError, "num_permutations"):
            IndexManifest.read(self.dir)

    def test_remove(self):
        IndexManifest.remove(self.dir)  # no manifest: no error
        self.manifest.write(self.dir)
        IndexManifest.remove(self.dir)
        self.assertFalse((self.dir / MANIFEST).exists())


if __name__ == "__main__":
    unittest.main()
