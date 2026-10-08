import unittest
from pathlib import Path
from tempfile import TemporaryDirectory

from hogprof.utils.phylo import TreeValidator


class TreeValidatorTest(unittest.TestCase):
    def setUp(self):
        self.temporary_directory = TemporaryDirectory()
        self.tmp_path = Path(self.temporary_directory.name)

    def tearDown(self):
        self.temporary_directory.cleanup()

    def write_tree(self, newick: str, suffix=".nwk") -> Path:
        path = self.tmp_path / f"species{suffix}"
        path.write_text(newick)
        return path

    def test_repairs_internal_species(self):
        tree_path = self.write_tree(
            "(((strain_a:1,strain_b:1)species_a:1,leaf_b:1)species_b:1,"
            "leaf_c:1)root;"
        )

        tree = TreeValidator(
            tree_path, {"species_a", "species_b", "leaf_c"}
        ).run()

        for name in ("species_a", "species_b", "leaf_c"):
            nodes = tree.search_nodes(name=name)
            self.assertEqual(1, len(nodes))
            self.assertTrue(nodes[0].is_leaf())

        self.assertEqual(
            {
                "strain_a",
                "strain_b",
                "species_a",
                "leaf_b",
                "species_b",
                "leaf_c",
            },
            {leaf.name for leaf in tree.iter_leaves()},
        )

    def test_rejects_missing_species(self):
        tree_path = self.write_tree("(species_a:1,species_b:1)root;")

        with self.assertRaisesRegex(ValueError, "missing species: species_c"):
            TreeValidator(
                tree_path, {"species_a", "species_b", "species_c"}
            ).run()

    def test_rejects_duplicate_species(self):
        tree_path = self.write_tree("(species_a:1,species_a:1)root;")

        with self.assertRaisesRegex(ValueError, "multiple tree nodes: species_a"):
            TreeValidator(tree_path, {"species_a"}).run()

    def test_rejects_non_newick(self):
        tree_path = self.write_tree(
            "(species_a:1,species_b:1)root;", suffix=".txt"
        )
        with self.assertRaisesRegex(ValueError, "newick format"):
            TreeValidator(tree_path, {"species_a", "species_b"}).run()


if __name__ == "__main__":
    unittest.main()
