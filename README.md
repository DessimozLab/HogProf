# HogProf
HogProf is an extensible and tunable approach to phylogenetic profiling using orthology data. It is powered by minhash-based data structures and computationally efficient.

# Features
  - Using orthoxoml files and a taxonomy calculated enhanced phylogenies of each family
  - These are transformed into minhash signatures and a locally sensitive hashing forest object for search and comparison of profiles
  - Taxonomic levels and evolutionary event types (presence, loss, duplication) can have custom weight in profile construction
  - Optimization of weights using machine learning

If you run into any problems feel free to contact me at [dmoi@unil.ch](dmoi@unil.ch)

# Quickstart
## Install from PyPI (recommended)

```bash
pip install hogprof
```

## Install from sources
```bash
git clone https://github.com/DessimozLab/HogProf.git
cd HogProf
pip install .
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
$python lshbuilder.py --outpath YourHogProfDirectory --dbtype all --OMA YourOmaDirectory/OmaServer.h5 --nthreads numberOfCPUcores         

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

Once the database is completed it can be interogated using a profiler object. Construction and usage of this object should be done using a python script or notebook. This shown in the example notebook `searchenrich.ipynb` found in the examples. Please feel free to modify it to suit the needs of your own research.
