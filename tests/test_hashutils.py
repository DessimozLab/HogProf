import unittest
from utils import hashutils


class HashutilsTest(unittest.TestCase):

    def test_hogid2fam(self):
        self.assertEqual(hashutils.hogid2fam('HOG:0279153'), 279153)

    def test_fam2hogid(self):
        self.assertEqual(hashutils.fam2hogid(279153), 'HOG:0279153')
        self.assertEqual(hashutils.fam2hogid('279153'), 'HOG:0279153')

    def test_result2hogid(self):
        self.assertEqual(hashutils.result2hogid('279153-gain'), 'HOG:0279153')
        self.assertEqual(hashutils.result2hogid('279153-gain-duplication'), 'HOG:0279153')

    def test_result2fam(self):
        self.assertEqual(hashutils.result2fam('279153-gain-duplication'), 279153)
        self.assertEqual(hashutils.result2fam('279153-gain'), 279153)

    def test_result2events(self):
        self.assertEqual(hashutils.result2events('279153-gain'), ['gain'])
        self.assertEqual(hashutils.result2events('279153-gain-duplication-loss'), ['gain', 'duplication', 'loss'])


if __name__ == "__main__":
    unittest.main()