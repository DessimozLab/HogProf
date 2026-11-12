"""
Defines the files produced by HogProf
"""

# contains Jaccard null approximations per level and thresholds
CALIBRATION = "calibration.json"

# LSH forest, answers "hog_key -> hog_keys with close signatures"
FOREST = "newlshforest.h5"

# Hash signature storage, answers "hog_key -> minhash signature"
# Required to compute Jaccard, as LSH doesn't store signatures
HASHSTORE = "hashes.h5"

# Main artifact of the output folder, defines the folder
# as HogProf output, contains metadata
MANIFEST = "manifest.json"

# Main metadata table, answers "hog_key -> all metadata of it"
MAPPING = "fam2orthoxml.csv"

MINHASH_GENERATOR = "wmg.pkl"
TAXA_INDEX = "taxaIndex.pkl"
TREE = "master_tree.corrected.nwk"
