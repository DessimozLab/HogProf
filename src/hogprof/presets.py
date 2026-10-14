"""Taxonomic ranges shared by the CLI and programmatic profile builder."""

DB_PRESETS = {
    "all": {"exclude_taxa": None, "taxon_mask": None},
    "eukarya": {"exclude_taxa": None, "taxon_mask": "Eukaryota"},

    # These below are the old masks  (pre 0.2.0) that use
    # the hardcoded taxa IDs. They should be substituted
    # with correct taxon names. To verify and enable later

    #'plants': {'exclude_taxa': None, 'taxon_mask': '33090'},
    #'archaea': {'exclude_taxa': None, 'taxon_mask': '2157'},
    #'bacteria': {'exclude_taxa': None, 'taxon_mask': '2'},
    #'protists': {
    #    'exclude_taxa': ['2', '2157', '33090', '4751', '33208'],
    #    'taxon_mask': None,
    # },
    #'fungi': {'exclude_taxa': None, 'taxon_mask': '4751'},
    #'metazoa': {'exclude_taxa': None, 'taxon_mask': '33208'},
    #'vertebrates': {'exclude_taxa': None, 'taxon_mask': '7742'},
}
