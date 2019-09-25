# HogProf

  - PyProfiler is an extensible and tunable approach to phylogenetic profiling using orthology data. It is powered by minhash based datastructures and computationally efficient.
  - Still under major development and may change
  - Magic

# Features

  - Using orthoxoml files and a taxonomy calculated enhanced phylogenies of each family
  - These are transformed into minhash signatures and a locally sensitive hashing forest object for search and comparison of profiles
  - Taxonomic levels and evolutionary event types ( presence, loss, duplication ) can have custom weight in profile construction
  - Optimization of weights using machine learning

If you run into any problems feel free to contact me at [dmoi@unil.ch](dmoi@unil.ch)

# Quickstart

```
$ git clone https://github.com/DessimozLab/pyprofiler.git
$ cd pyprofiler/pyprofiler
$ pip install -r req.txt .
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
$ mkdir YourPyProfilerDirectory
```
Now navigate to the pyprofiler source folder. Open the file config_utils.py in the utils folder and give it the location of you OMA data as well as the folder where you would like to save your pyprofiler databases.

```
$ cd utils
$ nano config_utils.py
```

change these to your parameters. Don't forget the trailing slash on your paths to your directories
```
config = {
    "dir":{
    "datadir": "YOURPYPROFILERDATADIRECTORY/",
    "omadir": "YOUROMADIRECTORY/"
    },
    "orthoxmltar":"",
    "email": "YOUREMAIL"
}
```
Your email will be used to identify you to the NCBI when using their API.

Ok. We're ready! Now let's compile a database containing all HOGs and our desired taxonomic levels using default settings. Launch the lshbuilder script from the pyprofiler folder.

dbtypes available on the command line are : all , plants , archaea, bacteria , eukarya , protists , fungi , metazoa and vertebrates.

```
$python lshbuilder.py --name YOURDBNAME --dbtype all                     
```

This should build a taxonomic tree for the genomes contained in the release and then calculate enhanced phylogenies for all HOGs in OMA.

Once the database is completed it can be interogated using a profiler object. Construction and usage of this object is shown in the example notebook searchenrich.ipynb found in the notebooks folder. It contains analysis related to a known and poorly described protein network. Please feel free to modify it to suit the needs of your own research.
