# HogProf

[![Tests](https://github.com/DessimozLab/HogProf/actions/workflows/ci.yml/badge.svg)](https://github.com/DessimozLab/HogProf/actions/workflows/ci.yml)
[![PyPI version](https://img.shields.io/pypi/v/HogProf.svg)](https://pypi.org/project/HogProf/)
[![DOI](https://img.shields.io/badge/DOI-10.1371%2Fjournal.pcbi.1007553-blue)](https://doi.org/10.1371/journal.pcbi.1007553)

HogProf is an extensible and tunable approach to phylogenetic profiling using orthology data. It is powered by minhash-based data structures and computationally efficient.

## Features

- Build phylogenetic profiles of families (HOGs) from orthology data.
- Profiles include ancestral states and evolutionary events, including gene losses and duplications.
- Search for similar profiles using locality-sensitive hashing (LSH).
- Generate profiles for subfamilies (subHOGs) to investigate evolutionary patterns at different taxonomic levels.

## Installation

### From PyPI (recommended)

```bash
pip install hogprof
```

### From source

```bash
git clone https://github.com/DessimozLab/HogProf.git
cd HogProf
pip install .
```

## Quickstart

HogProf provides two main commands:

- `hogprof build`: build a searchable index from an OMA database or OrthoXML files.
- `hogprof search`: search the index for HOGs with similar phylogenetic profiles.

To see all available options:

```bash
hogprof --help
hogprof build --help
hogprof search --help
```

## Building an index of phylogenetic profiles

HogProf can calculate phylogenetic profiles and build an index from either an OMA database or a collection of OrthoXML files. Both workflows require a species tree in Newick format, with node names matching the species identifiers in the input data.

### Option 1: Using an OMA database

Download the OMA HDF5 database from the [OMA Browser archives](https://omabrowser.org/oma/archives/).

**Note:** OMA database files can be very large. Check the size of the release before downloading.

For example, to download the July 2024 release:

```bash
mkdir oma_data
cd oma_data

wget https://omabrowser.org/All.Jul2024/OmaServer.h5
wget https://omabrowser.org/All.Jul2024/speciestree.nwk
```

Calculate the profiles and build the index:

```bash
hogprof build \
    --oma oma_data/OmaServer.h5 \
    --species-tree oma_data/speciestree.nwk \
    --output-dir hogprof_index \
    --jobs 20
```

Adjust `--jobs` to control the number of worker processes.

Use `--slice-subhogs` to generate profiles for subHOGs:

```bash
hogprof build \
    --oma oma_data/OmaServer.h5 \
    --species-tree oma_data/speciestree.nwk \
    --output-dir hogprof_index \
    --jobs 20 \
    --slice-subhogs 
```

### Option 2: Using FastOMA results

To build an index from FastOMA results, you need the FastOMA HOG OrthoXML file and the corresponding species tree.

**Step 1: Split the OrthoXML file**

You can do this using [`orthoxml-tools`](https://github.com/DessimozLab/orthoxml-tools). Install using `pip` in your preferred environment:

```bash
pip install orthoxml-tools
```

Activate the environment and split the FastOMA OrthoXML file into individual HOG files:

```bash
orthoxml-tools split \
    --infile /path/to/FastOMA_HOGs.orthoxml \
    --outdir /path/to/split
```

**Step 2: Build the HogProf index**

Provide the split OrthoXML files using `--orthoxml-glob` and the species tree used by FastOMA.

```bash
hogprof build \
    --orthoxml-glob "/path/to/split/*.orthoxml" \
    --species-tree "/path/to/species_tree_checked.nwk" \
    --output-dir hogprof_index \
    --jobs 20
```

Use `--slice-subhogs` to generate profiles for subHOGs:

```bash
hogprof build \
    --orthoxml-glob "/path/to/split/*.orthoxml" \
    --species-tree "/path/to/species_tree_checked.nwk" \
    --output-dir hogprof_index \
    --slice-subhogs \
    --jobs 20
```

Adjust `--jobs` to control the number of worker processes. Additional filtering options, such as `--min-species` and `--min-events`, are described under [Common build options](#common-build-options).

### Common build options

Run `hogprof build --help` for the complete list of options.

| Option | Description |
|---|---|
| `--oma` | Path to the OMA HDF5 database. |
| `--orthoxml-glob` | Glob expression matching OrthoXML files. Use this instead of `--oma`. |
| `--species-tree` | Newick species tree with node names matching the input data. |
| `--output-dir` | Directory where the index and other output files will be saved. |
| `--taxon-mask` | Restrict the analysis to a clade identified by its tree node name. |
| `--exclude-taxa` | Exclude specified clades. |
| `--min-species` | Minimum species threshold for retaining HOGs and subHOGs (default=10). |
| `--min-events` | Minimum number of loss or duplication events required to retain a profile (default=0). |
| `--db-type` | Taxonomic range preset: `all` or `eukarya`. |
| `--num-permutations` | Number of hash functions used to construct each profile (default=256). |
| `--seed` | Random seed for reproducibility (default=random). |
| `--loss-only` | Include only loss events. |
| `--duplication-only` | Include only duplication events. |
| `--use-tax-ids` | Match species to tree nodes using NCBI taxon IDs instead of species names. |
| `--slice-subhogs` | Generate profiles for subHOGs. |
| `--jobs` | Number of worker processes. |
| `--verbose` | Print additional progress information. |


## Searching an index

Once the index has been built, use it to search for HOGs with similar phylogenetic profiles.

For example:

```bash
hogprof search HOG:0000001 \
    --index-dir hogprof_index \
    --k 10
```

Replace `HOG:0000001` with the HOG identifier you want to query.

### Search options

| Option | Description |
|---|---|
| `hog_id` | HOG identifier to search for. |
| `--index-dir` | Directory containing the index generated by `hogprof build`. |
| `--k` | Number of nearest neighbors to consider. |
| `--level` | Taxonomic level at which the query is defined. Useful for ambiguous HOG identifiers occurring at multiple levels. |
| `--oma` | Path to the OMA HDF5 database. Necessary for indices build from OMA files. |
| `--nsamples` | Number of samples. Default=256. |
| `--slice-subhogs` | Enable subHOG search mode. |

For additional details, run:

```bash
hogprof search --help
```

## Additional analyses

HogProf also provides scripts and notebooks for downstream analyses, including:

- Tutorial showing basic functionalities.
- Preparing coevolution networks for visualization in Cytoscape.

These analyses require additional inputs beyond the index and search commands. See the relevant scripts and notebooks in the repository for details.

## Planned additions

- Generating a statistical threshold (empirical) to search for significantly similar profiles

## Citation

If you use HogProf in your research, please cite:

Moi D, Dessimoz C, Gonnet GH. (2020). *Efficient phylogenetic profiling using the OMA orthology database*. PLOS Computational Biology.

[https://doi.org/10.1371/journal.pcbi.1007553](https://doi.org/10.1371/journal.pcbi.1007553)

## Contact

For questions or issues, please open an issue in the [HogProf GitHub repository](https://github.com/DessimozLab/HogProf).