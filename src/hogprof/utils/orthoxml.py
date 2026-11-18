import collections
import xml.etree.ElementTree as ET
from io import StringIO


def strip_taxon_suffix(hog_id: str | None, taxon_id: str | None) -> str | None:
    """`HOG:F0783479.1a_32524` (taxonId 32524) -> `HOG:F0783479.1a`; other ids unchanged"""
    if hog_id and taxon_id and hog_id.endswith(f"_{taxon_id}"):
        return hog_id[:-len(taxon_id) - 1]
    return hog_id


class OrthoXMLBuilder:
    """
    Adds OMA-style subHOG IDs to OMA orthoxml (from Adrian): every
    orthologGroup gets `<hog name>_<taxonId>`, e.g. `HOG:F0783479.1a_32524`,
    and every geneRef the name of its HOG as `LOFT` (pyham names the levels
    it infers above a gene after it).

    Recent OMA releases ship these orthologGroup ids already: they are kept
    as they are, and the names are derived without their `_<taxonId>`.
    """

    def __init__(self):
        self.NS = "http://orthoXML.org/2011/"
        ET.register_namespace('', self.NS)  # Register default namespace

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
                # keep ids shipped by OMA, which use the same scheme
                if taxon_id is None:
                    # If no taxonId, just set the id to the og
                    if node.get('id') is None:
                        node.set('id', og)
                elif strip_taxon_suffix(node.get('id'), taxon_id) == node.get('id'):
                    node.set('id', og + "_" + taxon_id)
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
            # the root HOG's name, without the suffix if ids are already there
            og = strip_taxon_suffix(el.get('id'), el.get('taxonId')) or f"HOG:{i+1:08d}"
            dups = collections.defaultdict(int)
            rec_annotate(el, og)
        tree_str = StringIO()
        ### addition by Athina:
        tree_str.write('<?xml version="1.0" encoding="UTF-8"?>\n')
        doc.write(tree_str, encoding='unicode')
        return tree_str.getvalue()


def hog_names(orthoxml: str | ET.Element) -> dict[str, str]:
    """
    orthologGroup id -> HOG name, i.e. the id without its `_<taxonId>` suffix:
    `HOG:F0783479_117571` -> `HOG:F0783479`, `HOG:F0783479.1a_32524` -> `HOG:F0783479.1a`.
    The same name can be defined at several levels; (name, level) is unique.
    """
    xml = ET.fromstring(orthoxml) if isinstance(orthoxml, str) else orthoxml
    names = {}
    for group in xml.iter('{http://orthoXML.org/2011/}orthologGroup'):
        hog_id = group.get('id')
        if hog_id is not None:
            names[hog_id] = strip_taxon_suffix(hog_id, group.get('taxonId'))
    return names


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
