"""OMA HOG ids and names."""

import unittest
import xml.etree.ElementTree as ET

from hogprof.lshbuilder import LSHBuilder
from hogprof.utils import orthoxml as oxml

from split_family import TAXA, orthoxml

NS = "{http://orthoXML.org/2011/}"


def group_ids(xml):
    return [g.get("id") for g in ET.fromstring(xml).iter(NS + "orthologGroup")]


def loft_ids(xml):
    return {g.get("id"): g.get("LOFT") for g in ET.fromstring(xml).iter(NS + "geneRef")}


class OmaIdsTest(unittest.TestCase):
    def test_strip_taxon_suffix(self):
        self.assertEqual("HOG:F0783479.1a", oxml.strip_taxon_suffix("HOG:F0783479.1a_32524", "32524"))
        self.assertEqual("HOG:F0783479", oxml.strip_taxon_suffix("HOG:F0783479_117571", "117571"))
        # only its own taxonId is stripped
        self.assertEqual("HOG:F0783479_117571", oxml.strip_taxon_suffix("HOG:F0783479_117571", "32524"))
        self.assertEqual("HOG:F1", oxml.strip_taxon_suffix("HOG:F1", None))
        self.assertIsNone(oxml.strip_taxon_suffix(None, "1"))

    def test_shipped_ids_are_kept(self):
        shipped = orthoxml(with_ids=True)
        annotated = oxml.OrthoXMLBuilder().add_loft_ids(shipped)
        # previously the root became HOG:F1_1_1 and subHOGs HOG:F1_1.1a_3
        self.assertEqual(group_ids(shipped), group_ids(annotated))
        self.assertIn(f"HOG:F1_{TAXA['R']}", group_ids(annotated))

    def test_missing_ids_get_the_same_scheme(self):
        # OMA releases without ids: the same ids as OMA ships are generated
        generated = oxml.OrthoXMLBuilder().add_loft_ids(orthoxml(with_ids=False))
        self.assertEqual(group_ids(orthoxml(with_ids=True)), group_ids(generated))
        self.assertEqual(loft_ids(oxml.OrthoXMLBuilder().add_loft_ids(orthoxml(with_ids=True))),
                         loft_ids(generated))

    def test_genes_get_their_hog_name_as_loft(self):
        lofts = loft_ids(oxml.OrthoXMLBuilder().add_loft_ids(orthoxml()))
        self.assertEqual("HOG:F1.1a", lofts["0"])    # m1, first copy
        self.assertEqual("HOG:F1.2b", lofts["21"])   # d1, second copy
        self.assertEqual("HOG:F1", lofts["40"])      # y, not duplicated

    def test_hog_names(self):
        names = oxml.hog_names(oxml.OrthoXMLBuilder().add_loft_ids(orthoxml()))
        self.assertEqual("HOG:F1", names[f"HOG:F1_{TAXA['R']}"])
        self.assertEqual("HOG:F1.2b", names[f"HOG:F1.2b_{TAXA['D']}"])
        self.assertEqual({"HOG:F1", "HOG:F1.1a", "HOG:F1.1b", "HOG:F1.2a", "HOG:F1.2b"}, set(names.values()))

    def test_subhog_name(self):
        names = {"HOG:F1.1a_3": "HOG:F1.1a"}
        self.assertEqual("HOG:F1.1a", LSHBuilder._subhog_name("HOG:F1.1a_3", names))
        self.assertEqual("HOG:F1.1a__2", LSHBuilder._subhog_name("HOG:F1.1a_3__2", names))
        # FastOMA (orthoxml input): ids as they are
        self.assertEqual("HOG:0055821_3", LSHBuilder._subhog_name("HOG:0055821_3", {}))


if __name__ == "__main__":
    unittest.main()
