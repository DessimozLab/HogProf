"""Event-count filtering of tree profiles."""

import unittest

import ete3

from hogprof.utils.pyhamutils import _has_min_events, count_events


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


if __name__ == "__main__":
    unittest.main()
