import collections
import xml.etree.ElementTree as ET
from io import StringIO


class OrthoXMLBuilder:
    """
    This is a function for adding subHOG IDs to the OMA orthoxml files (from Adrian).
    Will become redundant with the next version of OMA
    """

    def __init__(self):
        self.NS = "http://orthoXML.org/2011/"
        ET.register_namespace('', self.NS)  # Register default namespace

    def add_loft_ids(self, orthoxml):

        def encode_paralog_cluster_id(prefix, nr):
            letters = []
            while nr // 26 > 0:
                letters.append(chr(97 + (nr % 26)))
                nr = nr // 26 - 1
            letters.append(chr(97 + (nr % 26)))
            return prefix + ''.join(letters[::-1])

        def next_sub_hog_id(idx):
            dups[idx] += 1
            return dups[idx]

        def rec_annotate(node, og, idx=0):
            if node.tag == f"{{{self.NS}}}orthologGroup":
                taxon_id = node.get('taxonId')
                if taxon_id is not None:
                    node.set('id', og + "_" + taxon_id)
                else:
                    # If no taxonId, just set the id to the og
                    node.set('id', og )
                for child in list(node):
                    rec_annotate(child, og, idx)
            elif node.tag == f"{{{self.NS}}}paralogGroup":
                idx += 1
                next_og = f"{og}.{next_sub_hog_id(idx)}"
                for i, child in enumerate(list(node)):
                    rec_annotate(child, encode_paralog_cluster_id(next_og, i), idx)
            elif node.tag == f"{{{self.NS}}}geneRef":
                # Set the id to the geneRef id
                node.set('LOFT', og)

        doc = ET.parse(StringIO(orthoxml))
        xml = doc.getroot()
        for i, el in enumerate(xml.findall(".//groups/orthologGroup", {"": self.NS})):
            og = el.get('id', f"HOG:{i+1:08d}")
            dups = collections.defaultdict(int)
            rec_annotate(el, og)
        tree_str = StringIO()
        ### addition by Athina:
        tree_str.write('<?xml version="1.0" encoding="UTF-8"?>\n')
        doc.write(tree_str, encoding='unicode')
        return tree_str.getvalue()


def count_species(ortho_fam):
    return ortho_fam.count('<species name=')


def iter_species(filename):
    root = ET.parse(filename).getroot()
    for child in root:
        if child.tag.rsplit("}", 1)[-1] == "species":
            yield child