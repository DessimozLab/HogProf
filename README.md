# HogProf
  - HogProf is an extensible and tunable approach to phylogenetic profiling using orthology data. It is powered by minhash based datastructures and computationally efficient.
  - Still under major development and may change

# Features

  - Using orthoxoml files and a taxonomy calculated enhanced phylogenies of each family
  - These are transformed into minhash signatures and a locally sensitive hashing forest object for search and comparison of profiles
  - Taxonomic levels and evolutionary event types ( presence, loss, duplication ) can have custom weight in profile construction
  - Optimization of weights using machine learning

If you run into any problems feel free to contact me at [dmoi@unil.ch](dmoi@unil.ch)

# Quickstart

to install from github
```
$ git clone https://github.com/DessimozLab/HogProf.git
$ pip install -r pipreqs.txt .
```

or install it from pypi with pip

```
$ pip install hogprof
```

lets get a current version of the OMA hdf5 file and GAF. This will alow us to use the HOGs and study the functional enrichment of our search results.


```
$ cd ../..
$ mkdir YourOmaDirectory
$ cd YourOmaDirectory
$ wget https://omabrowser.org/All/OmaServer.h5
$ wget https://omabrowser.org/All/oma-go.txt.gz
```

We also need to make a location to store our pyprofiler databases

```
$ cd ..
$ mkdir YourHogProfDirectory
```

Ok. We're ready! Now let's compile a database containing all HOGs and our desired taxonomic levels using default settings. Launch the lshbuilder.
dbtypes available on the command line are : all , plants , archaea, bacteria , eukarya , protists , fungi , metazoa and vertebrates. These will use the NCBI taxonomy as a tree to annotate events in different gene family's histories.
```
$python lshbuilder.py --outpath YourHogProfDirectory --dbtype all --OMA YourOmaDirectory/OmaServer.h5 --nthreads numberOfCPUcores          

```
This should build a taxonomic tree for the genomes contained in the release and then calculate enhanced phylogenies for all HOGs in OMA.

Once the database is completed it can be interogated using a profiler object. Construction and usage of this object should be done using a python script or notebook. This shown in the example notebook searchenrich.ipynb found in the examples. Please feel free to modify it to suit the needs of your own research.