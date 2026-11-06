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

    @staticmethod
    def _has_level_id(group: ET.Element) -> bool:
        """True if the group id is already in the '<HOG id>_<taxonId>' form"""
        group_id, taxon_id = group.get('id'), group.get('taxonId')
        if not group_id:
            return False
        return taxon_id is None or group_id.endswith(f"_{taxon_id}")

    def add_loft_ids(self, orthoxml: str | ET.Element) -> str:
        """Annotate XML text or an already parsed root, modifying it in place."""

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

        xml = ET.fromstring(orthoxml) if isinstance(orthoxml, str) else orthoxml
        doc = ET.ElementTree(xml)
        for i, el in enumerate(xml.findall(".//groups/orthologGroup", {"": self.NS})):
            # Newer OMA releases (e.g. All.May2026) already give every subHOG
            # its id, in the same "<HOG id>_<taxonId>" form. Keep those:
            # annotating again would turn them into ids OMA does not know
            if all(self._has_level_id(group) for group in el.iter(f"{{{self.NS}}}orthologGroup")):
                continue
            og = el.get('id', f"HOG:{i+1:08d}")
            dups = collections.defaultdict(int)
            rec_annotate(el, og)
        tree_str = StringIO()
        ### addition by Athina:
        tree_str.write('<?xml version="1.0" encoding="UTF-8"?>\n')
        doc.write(tree_str, encoding='unicode')
        return tree_str.getvalue()


def get_root_hog_id(orthoxml: str | ET.Element) -> str:
    xml = ET.fromstring(orthoxml) if isinstance(orthoxml, str) else orthoxml
    groups = xml.findall(
        '{http://orthoXML.org/2011/}groups/{http://orthoXML.org/2011/}orthologGroup'
    )
    if len(groups) != 1:
        raise ValueError(f'Expected one root HOG, found {len(groups)}')
    return groups[0].attrib['id']


def count_species(ortho_fam: str | ET.Element) -> int:
    if isinstance(ortho_fam, ET.Element):
        return len(ortho_fam.findall('{http://orthoXML.org/2011/}species'))
    return ortho_fam.count('<species name=')


def iter_species(filename):
    root = ET.parse(filename).getroot()
    for child in root:
        if child.tag.rsplit("}", 1)[-1] == "species":
            yield child
