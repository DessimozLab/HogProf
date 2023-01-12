import unittest
from utils import files_utils


class HpputilsTest(unittest.TestCase):

    def test_replace_characters(self):
        self.assertEqual(files_utils.replace_characters('test.replace,character(to underscore)working:well'),
                         'test_replace_character_to_underscore_working_well')


if __name__ == "__main__":
    unittest.main()