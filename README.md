# HogProf
HogProf is an extensible and tunable approach to phylogenetic profiling using orthology data. It is powered by minhash-based data structures and computationally efficient.

⚠️ Development branch

The implementation supporting taxonomic-level-aware phylogenetic profiling is currently available on the levels_addition branch and has not yet been merged into master.

If you would like to use the latest functionality described here, make sure you clone or install the levels_addition branch rather than the default master branch.

# Features
  - Using orthoxoml files and a taxonomy calculated enhanced phylogenies of each family
  - These are transformed into minhash signatures and a locally sensitive hashing forest object for search and comparison of profiles
  - Taxonomic levels and evolutionary event types (presence, loss, duplication) can have custom weight in profile construction
  - Optimization of weights using machine learning

If you run into any problems, feel free to contact me at [dmoi@unil.ch](dmoi@unil.ch)

# Quickstart
## Install from PyPI (recommended)

<<<<<<< HEAD
to install from github using a conda environment (recommended)
||||||| 01f5de5
to install from github
=======
```bash
pip install hogprof
>>>>>>> develop
```
<<<<<<< HEAD
$ git clone --branch levels_addition https://github.com/DessimozLab/HogProf.git
$ cd HogProf
$ conda create -n hogprof_levels python=3.8 pip
$ conda activate hogprof_levels
$ pip install -r pipreqs.txt
$ pip install -e .
```
or to install from pypi (levels mode not supported)
```
$ pip install hogprof
||||||| 01f5de5
$ git clone https://github.com/DessimozLab/HogProf.git
$ pip install -r pipreqs.txt .
```
or to install from pypi
```
$ pip install hogprof
=======

## Install from sources
```bash
git clone https://github.com/DessimozLab/HogProf.git
cd HogProf
pip install .
>>>>>>> develop
```


# Usage

## Example: using the OMA database
Let's get a current version of the OMA hdf5 file and GAF. This will aloww us to use the HOGs and study the functional enrichment of our search results. **Careful, this download is heavy (hundreds of GB)**:

```bash
mkdir YourOmaDirectory
cd YourOmaDirectory
wget https://omabrowser.org/All/OmaServer.h5
wget https://omabrowser.org/All/oma-go.txt.gz
```

The latest release (May.2026) is more than 300GB in size. Alternatively, you can use earlier releases, for example:
```bash
wget https://omabrowser.org/All.Jul2024/OmaServer.h5
```
which is 180GB in size. For other releases, check `https://omabrowser.org/oma/archives/`, select release and search for "OMA Browser database (as hdf5)". 


We also need to make a location to store our pyprofiler databases

```bash
cd ..
mkdir YourHogProfDirectory
```

Let's now compile a database containing all HOGs and our desired taxonomic levels using default settings. Launch the lshbuilder as shown below.
`dbtypes` available on the command line are: `all`, `plants`, `archaea`, `bacteria`, `eukarya`, `protists`, `fungi`, `metazoa`, and `vertebrates`. These will use the NCBI taxonomy as a tree to annotate events in different gene family's histories.

If you are using an OMA release before 2022 you will need to use the NCBI tree. This is the default tree used by HogProf.

```
lshbuilder --outpath YourHogProfDirectory --dbtype all --OMA YourOmaDirectory/OmaServer.h5 --nthreads numberOfCPUcores         

```

If you are using OMA releases after 2022 we will also be downloading the OMA taxonomic tree since it is more accurate than the NCBI tree. This will be used to annotate the events in the gene family histories.

```bash
wget https://omabrowser.org/All/speciestree.nwk
```

We are ready to use the OMA tree to build our database.

```bash
lshbuilder --outpath YourHogProfDirectory --dbtype all --OMA YourOmaDirectory/OmaServer.h5 --nthreads numberOfCPUcores --mastertree YourOmaDirectory/speciestree.nwk --reformat_names True

```

This should build a taxonomic tree for the genomes contained in the release and then calculate enhanced phylogenies for all HOGs in OMA.

<<<<<<< HEAD
Once the database is completed it can be interogated using a profiler object. Construction and usage of this object should be done using a python script or notebook. This shown in the example notebook searchenrich.ipynb found in the examples. Please feel free to modify it to suit the needs of your own research.


# levels mode

The installation process is the same (github only for now).

A fastOMA result folder is expected as input.

**Step 1**: split the orthoxml into HOG-specific orthoxmls, using this script (requires a fastOMA installation):
``` 
python path/to/HogProf/src/HogProf/splitter.py \
    /path/to/output_directory \
    --fastoma-path /path/to/FastOMA
```

**Step 2**: build the HogProf database in **levels mode** by running `lshbuilder.py` with the `--slicesubhogs` option enabled. Point `--OrthoGlob` to the HOG-specific OrthoXML files generated in Step 1 and provide the species tree used by FastOMA.

```bash
python path/to/HogProf/src/HogProf/lshbuilder.py \
    --OrthoGlob "/path/to/fastoma_results/splits/*.orthoxml" \
    --outpath "/path/to/output_database/" \
    --mastertree "/path/to/fastoma_results/species_tree_checked.nwk" \
    --reformat_names True \
    --slicesubhogs True \
    --eventslim 1 \
    --nthreads 20
```

Required parameters:
- `--OrthoGlob`: glob matching the HOG-specific OrthoXML files created in Step 1.
- `--outpath`: output directory for the HogProf database.
- `--mastertree`: species tree corresponding to the FastOMA run.
- `--slicesubhogs True`: enables levels mode by generating profiles for subHOGs.

Common optional parameters:
- `--eventslim`: minimum number of evolutionary events required to retain a subHOG (e.g. `1`).
- `--specieslim`: minimum number of species required in a subHOG (default: `10`).
- `--nthreads`: number of CPU threads to use.
- `--reformat_names True`: recommended if the species tree labels need to be converted to numeric identifiers.

**Step 3**: generate Jaccard similarity distributions and compute threshold tables from the HogProf database.

This step analyzes the output of `lshbuilder` and produces the **bin-wise Jaccard threshold table** used for downstream filtering and scoring.

```bash
python path/to/HogProf/src/HogProf/generate_jaccard_distr.py \
    --hogprof_folder "/path/to/2a_hogprof_output/" \
    --output "/path/to/output_distributions/" \
    --profiler-path "/path/to/profiler.py"
```

Required parameters:
- `--hogprof_folder`: folder containing the HogProf LSHbuilder output (Step 2).
- `--output`: directory where distribution plots and threshold tables will be written.

Optional parameters:
- `--profiler-path`: path to `profiler.py`.  
  By default, this points to the `profiler.py` located in the same directory as the script.

### Output

This step produces a key file used for thresholding:

- `bins_jaccard_thresholds.csv`  
  → table of Jaccard similarity thresholds computed per bin, used for downstream scoring and filtering.

The file is saved in the directory specified by `--output`.

**Step 4**: extract level-specific hits using the HogProf profiler.

This step runs `extract_levels_hits.py` to query the HogProf database built in Step 2 and generate a ranked list of hits per subHOG using the threshold table generated in Step 3.

```bash
python path/to/HogProf/src/HogProf/extract_levels_hits.py \
    --input "/path/to/hogprof_output/" \
    --outputfile "/path/to/output/extracted_hits.csv" \
    --queries_file "/path/to/queries.tsv" \
    --thresholds "/path/to/bins_jaccard_thresholds.csv" \
    --k 1000 \
    --allvsall True \
    --profilerpath "profiler.py"
```

Required parameters:
- `--input`: HogProf LSHbuilder output directory (Step 2).  
  Must contain files such as:
  `newlshforest.pkl`, `hashes.h5`, `reformatted_tree.nwk`, `fam2orthoxml.csv`, `profilersavingpath.csv`.
- `--outputfile`: path to the CSV file where extracted hits will be saved.
- `--queries_file`: table of query subHOGs (e.g. venom HOGs).  
  Must contain a `hashid` column and a subHOG identifier column.
- `--thresholds`: CSV file generated in Step 3 (`bins_jaccard_thresholds.csv`).

Optional parameters:
- `--k`: number of top hits to retain per query (default: `1000`).
- `--allvsall`: also performs all-vs-all comparisons (default: `True`, may be computationally expensive).
- `--profilerpath`: path to `profiler.py` from HogProf (default assumes local availability).

### Output

This step produces:

- `extracted_hits.csv`  
  → ranked hit table containing similarity matches between query subHOGs and the HogProf database, filtered using bin-specific Jaccard thresholds.

**Step 5**: visualize coevolution networks in Cytoscape.

This final step uses the Jupyter notebook `cytoscape_import.ipynb` to generate and export coevolution networks from the extracted HogProf hits.

The notebook loads the output from Step 4 and prepares network files that can be visualized in Cytoscape.

In addition to the extracted hits, the notebook requires a **nodes metadata table** corresponding to the query subHOGs used in Step 4.

||||||| 01f5de5
Once the database is completed it can be interogated using a profiler object. Construction and usage of this object should be done using a python script or notebook. This shown in the example notebook searchenrich.ipynb found in the examples. Please feel free to modify it to suit the needs of your own research.
=======
Once the database is completed it can be interogated using a profiler object. Construction and usage of this object should be done using a python script or notebook. This shown in the example notebook `searchenrich.ipynb` found in the examples. Please feel free to modify it to suit the needs of your own research.
>>>>>>> develop
