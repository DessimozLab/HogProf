HogProf
=====================

  - HogProf is an extensible and tunable approach to phylogenetic profiling using orthology data. It is powered by minhash based datastructures and computationally efficient.
  - Still under major development and may change.

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   installation
   usage
   troubleshooting
   credits

Installation
------------

To install My Project, run the following command:

Using pip
.. code-block:: bash
   $ pip install hogprof

Or from github 
.. code-block:: bash
$ git clone https://github.com/DessimozLab/HogProf.git
$ pip install -r pipreqs.txt .

Quickstart
-----

To use the library we need an OMA instance's HDF5 file containing HOG info and some accesory files.

.. code-block:: bash
   $ mkdir YourOmaDirectory
   $ cd YourOmaDirectory
   $ wget https://omabrowser.org/All/OmaServer.h5
   $ wget https://omabrowser.org/All/oma-go.txt.gz

Let's create a directory for the phylogenetic rpfiling database were going to make.


.. code-block:: bash
   $ mkdir YourDBDirectory

Ok. We're ready! Now let's compile a database containing all HOGs and our desired taxonomic levels using default settings. Launch the lshbuilder.
dbtypes available on the command line are : all , plants , archaea, bacteria , eukarya , protists , fungi , metazoa and vertebrates. These will use the NCBI taxonomy as a tree to annotate events in different gene family's histories.

.. code-block:: bash
   $python lshbuilder.py --outpath YourHogProfDirectory --dbtype all --OMA YourOmaDirectory/OmaServer.h5 --nthreads numberOfCPUcores          

This should build a taxonomic tree for the genomes contained in the release and then calculate enhanced phylogenies for all HOGs in OMA.
Once the database is completed it can be interogated using a profiler object. Construction and usage of this object should be done using a python script or notebook.


.. code-block:: python
   import HogProf
   myproject.do_x()




Troubleshooting
---------------


If you encounter any issues while using My Project, please file a bug report on our GitHub repository: https://github.com/user/repo/issues


Credits
-------

My Project was created by John Doe.