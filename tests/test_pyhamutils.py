"""Event-count filtering of tree profiles."""

import unittest

import ete3
import pyham

from hogprof.utils import orthoxml as oxml
from hogprof.utils.pyhamutils import _has_min_events, count_events, get_ham_treemap_from_row, merge_split_hogs

from split_family import SPECIES, TREE, orthoxml


def treemap(events):
    """A tree profile with (duplications, losses) per node"""
    tree = ete3.Tree("((a,b)ab,c)root;", format=1)
    events = list(events)
    for i, node in enumerate(tree.traverse()):
        dupl, lost = events[i] if i < len(events) else (0, 0)
        node.add_features(dupl=dupl, lost=lost)
    return tree


class MinEventsTest(unittest.TestCase):
    def test_counts_losses_and_duplications(self):
        self.assertEqual(4, count_events(treemap([(1, 0), (0, 2), (1, 0)])))

    def test_zero_keeps_profiles_without_events(self):
        self.assertTrue(_has_min_events(treemap([]), 0))
        self.assertTrue(_has_min_events(treemap([(0, 0)] * 5), 0))

    def test_minimum_is_inclusive_over_all_events(self):
        # one duplication + one loss = two events, neither type alone reaches 2
        profile = treemap([(1, 0), (0, 1)])
        self.assertTrue(_has_min_events(profile, 2))
        self.assertFalse(_has_min_events(profile, 3))


class SplitHogsTest(unittest.TestCase):
    """pyham splits the gene at C into one HOG per duplication; we merge them"""

    def ham(self):
        xml = oxml.OrthoXMLBuilder().add_loft_ids(orthoxml())
        ham = pyham.Ham(TREE, xml, type_hog_file="orthoxml", tree_format="newick_string",
                        use_internal_name=True, orthoXML_as_string=True)
        return ham, ham.get_list_top_level_hogs()[0]

    @staticmethod
    def hogs_at(root, level):
        return [h for h in root.get_all_descendant_hogs() if h.genome.name == level]

    def test_pyham_splits(self):
        # documents the pyham behaviour the merge exists for
        _, root = self.ham()
        split = self.hogs_at(root, "C")
        self.assertEqual(2, len(split))
        self.assertEqual({4}, {len(h.get_all_descendant_genes()) for h in split})
        self.assertEqual(1, len({h.hog_id for h in split}))

    def test_merge(self):
        ham, root = self.ham()
        self.assertEqual(1, merge_split_hogs(root))

        merged, = self.hogs_at(root, "C")
        self.assertEqual(8, len(merged.get_all_descendant_genes()))
        self.assertEqual(2, len(merged.duplications))
        self.assertTrue(all(d.parent is merged for d in merged.duplications))
        self.assertEqual([merged], [g for g in merged.genome.genes if g.parent is root])
        # copies keep their duplication
        self.assertEqual(4, sum(bool(c.arose_by_duplication) for c in merged.children))

        # nothing left to merge
        self.assertEqual(0, merge_split_hogs(root))

    def test_profiles_have_no_spurious_losses(self):
        row = (0, oxml.OrthoXMLBuilder().add_loft_ids(orthoxml()))
        kwargs = dict(swap_ids=False, orthoXML_as_string=True, limit_species=1)

        levels = get_ham_treemap_from_row(row, TREE, slicesubhogs=True, **kwargs)
        self.assertFalse([key for key in levels if "__" in key])
        at_c, = [tm for key, tm in levels.items() if key.startswith("C_")]

        for profile in (at_c, get_ham_treemap_from_row(row, TREE, **kwargs)):
            nodes = {n.name: n for n in profile.traverse()}
            self.assertEqual(0, sum(n.lost or 0 for n in profile.traverse()))
            self.assertEqual(1, nodes["C"].nbr_genes)
            self.assertEqual((2, 2), (nodes["M"].nbr_genes, nodes["D"].nbr_genes))
            self.assertEqual(sorted(SPECIES[:4]), sorted(n.name for n in at_c.get_leaves()))


if __name__ == "__main__":
    unittest.main()
