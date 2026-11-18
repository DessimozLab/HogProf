"""
A small family where pyham splits a HOG, shaped like HOG:F0412788 of OMA May2026.

Species tree: R -> (C -> (M -> m1, m2), (D -> d1, d2)), y
The root HOG is at R. Two separate duplications happened between C and the
clades M and D, so pyham places both at C, but the orthoxml has no group at C
(it only lists levels where something changes): pyham then inserts one HOG at
C per duplication, splitting the single ancestral gene at C in two.

Group ids follow OMA (`<name>_<taxonId>`), as shipped since OMA May2026.
"""

TREE = "((((m1,m2)M,(d1,d2)D)C,y)R);"

TAXA = {"R": 1, "C": 2, "M": 3, "D": 4}

SPECIES = ["m1", "m2", "d1", "d2", "y"]


def _species():
    lines = []
    for i, name in enumerate(SPECIES):
        genes = "".join(f'<gene id="{10 * i + j}" protId="{name}_{j}"/>' for j in range(4))
        lines.append(
            f'<species name="{name}" NCBITaxId="{100 + i}">'
            f'<database name="db" version="1"><genes>{genes}</genes></database></species>'
        )
    return "\n".join(lines)


def _group(name, level, genes, with_ids, root=False):
    attributes = f' taxonId="{TAXA[level]}"'
    if with_ids:
        attributes += f' id="{name}_{TAXA[level]}"'
    elif root:
        # before May2026, OMA gave only the root a (bare) id
        attributes += f' id="{name}"'
    refs = "".join(f'<geneRef id="{g}"/>' for g in genes)
    return f'<orthologGroup{attributes}><property name="TaxRange" value="{level}"/>{refs}</orthologGroup>'


def orthoxml(with_ids=True):
    """The family; with_ids=False: no ids but the bare root id, like OMA releases before May2026"""
    # gene ids: species index * 10 + copy
    m1, m2, d1, d2, y = (10 * i for i in range(5))
    groups = (
        _group("HOG:F1", "R", [], with_ids, root=True).replace(
            "</orthologGroup>",
            "<paralogGroup>"
            + _group("HOG:F1.1a", "M", [m1, m2], with_ids)
            + _group("HOG:F1.1b", "M", [m1 + 1, m2 + 1], with_ids)
            + "</paralogGroup>"
            + "<paralogGroup>"
            + _group("HOG:F1.2a", "D", [d1, d2], with_ids)
            + _group("HOG:F1.2b", "D", [d1 + 1, d2 + 1], with_ids)
            + "</paralogGroup>"
            + f'<geneRef id="{y}"/>'
            + "</orthologGroup>",
        )
    )
    return (
        '<?xml version="1.0" encoding="UTF-8"?>\n'
        '<orthoXML xmlns="http://orthoXML.org/2011/" origin="test" originVersion="1" version="0.5">\n'
        f"{_species()}\n<groups>{groups}</groups></orthoXML>"
    )
