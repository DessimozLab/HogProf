import argparse
import functools
import glob
import pandas as pd
import time as t
import pickle
import xml.etree.cElementTree as ET
from ete3 import Phyloxml
import sys
import traceback
from datasketch import MinHashLSHForest , WeightedMinHashGenerator
from datetime import datetime
import h5py
import time
import gc
from pyoma.browser import db
from HogProf.utils import pyhamutils, hashutils , files_utils
import numpy as np
import tqdm
import random
import tqdm
import logging
import multiprocessing as mp
import os
import ete3
import collections
import io
import pickle
import queue
import random
import sys
import time
import time as t
import xml.etree.ElementTree as ET
from datetime import datetime
from pathlib import Path
from HogProf import __version__

import h5py
import numpy as np
import pandas as pd
import tqdm
from datasketch import MinHashLSHForest, WeightedMinHashGenerator
from pyoma.browser import db
from tables import open_file

from HogProf.utils import hashutils, phylo, pyhamutils

logger = logging.getLogger(__name__)

random.seed(0)
np.random.seed(0)
#import psutil
#process = psutil.Process()

'''This is a function for adding subhog ids to the OMA orthoxml files (from Adrian). Will become redundant
 with the next version of OMA'''
class OrthoXMLBuilder:

    def __init__(self):
        self.NS = "http://orthoXML.org/2011/"
        ET.register_namespace('', self.NS)  # Register default namespace

    def add_loft_ids(self, orthoxml):

        def encodeParalogClusterId(prefix, nr):
            letters = []
            while nr // 26 > 0:
                letters.append(chr(97 + (nr % 26)))
                nr = nr // 26 - 1
            letters.append(chr(97 + (nr % 26)))
            return prefix + ''.join(letters[::-1])

        def nextSubHogId(idx):
            dups[idx] += 1
            return dups[idx]

        def rec_annotate(node, og, idx=0):
            if node.tag == f"{{{self.NS}}}orthologGroup":
                taxon_id = node.get('taxonId')
                if taxon_id is not None:
                    node.set('id', og + "_" + taxon_id)
                else:
                    # If no taxonId, just set the id to the og
                    node.set('id', og )
                for child in list(node):
                    rec_annotate(child, og, idx)
            elif node.tag == f"{{{self.NS}}}paralogGroup":
                idx += 1
                next_og = "{}.{}".format(og, nextSubHogId(idx))
                for i, child in enumerate(list(node)):
                    rec_annotate(child, encodeParalogClusterId(next_og, i), idx)
            elif node.tag == f"{{{self.NS}}}geneRef":
                # Set the id to the geneRef id
                node.set('LOFT', og)

        doc = ET.parse(io.StringIO(orthoxml))
        xml = doc.getroot()
        for i, el in enumerate(xml.findall(".//groups/orthologGroup", {"": self.NS})):
            og = el.get('id', "HOG:{:08d}".format(i+1))
            dups = collections.defaultdict(int)
            rec_annotate(el, og)
        tree_str = io.StringIO()
        ### addition by Athina:
        tree_str.write('<?xml version="1.0" encoding="UTF-8"?>\n')
        doc.write(tree_str, encoding='unicode')
        return tree_str.getvalue()



class LSHBuilder:
    """
    This class contains the stuff you need to make 
    a phylogenetic profiling 
    database with input orthxml files and a taxonomic tree
    You must either input an OMA hdf5 file or an ensembl tarfile 
    containing orthoxml file with orthologous groups.

    You can provide a species tree or use the ncbi taxonomy 
    with a list of taxonomic codes for all the species in your db
    """

    def __init__(self,h5_oma=None,fileglob = None, taxa=None,masterTree=None, saving_name=None ,   numperm = 256,  treeweights= None , taxfilter = None, taxmask= None , lossonly = False, duplonly = False, verbose = False , use_taxcodes = False , datetime = datetime.now() , reformat_names = False,
                  slicesubhogs=False, limit_species=10, limit_events=0):
                
        """
            Initializes the LSHBuilder class with the specified parameters and sets up the necessary objects.
            
            Args:
            - tarfile_ortho (str):  path to an ensembl tarfile containing orthoxml files
            - h5_oma (str): path to an OMA hdf5 file
            - masterTree (str): path to a newick tree file
            - saving_name (str): path to the directory where the output files will be saved
            - numperm (int): the number of permutations to use in the MinHash generation (default: 256)
            - treeweights (str): path to a pickled file containing the weights for the tree
            - taxfilter (str): path to a file containing a list of taxonomic codes to filter from the tree
            - taxmask (str): path to a file containing a list of taxonomic codes to mask from the tree
            - verbose (bool): whether to print verbose output (default: False)
            - slicesubhogs (bool): whether to slice the subhogs (default: False)
            - limit_species (int): the minimum number of species in a subHOG that is included in the database (default: 10)
            - limit_events (int): the minimum number of events (loss,duplication) in a subHOG that is included in the database (default: 0)

        """
        print("\nInitializing LSHBuilder")

        if h5_oma:
            self.h5OMA = h5_oma
            self.db_obj = db.Database(h5_oma)
            self.oma_id_obj = db.OmaIdMapper(self.db_obj)
        else:
            self.h5OMA = None
            self.db_obj = None
            self.oma_id_obj = None
        
        self.reformat_names = reformat_names
        self.slicesubhogs = slicesubhogs
        self.swap2taxcode = use_taxcodes
        self.tax_filter = taxfilter
        self.tax_mask = taxmask
        self.verbose = verbose
        self.datetime = datetime
        self.fileglob = fileglob
        self.idmapper = None
        self.date_string = "{:%B_%d_%Y_%H_%M}".format(datetime.now())
        self.limit_species = limit_species
        self.limit_events = limit_events
        if saving_name:
            self.saving_name= saving_name 
            if self.saving_name[-1]!= '/':
                self.saving_name = self.saving_name+'/'
            self.saving_path = saving_name
            if not os.path.isdir(self.saving_path):
                os.mkdir(path=self.saving_path)
        else:
            raise Exception( 'please specify an output location' )

        species = self._get_species_names()
        if masterTree is None:
            if not h5_oma:
                raise TypeError('Please specify either a database or a tree')

            self.tree_string, self.tree = phylo.get_tree(genomes=species, outdir=self.saving_path)
        else:
            # validate tree
            self.tree = phylo.TreeValidator(masterTree, species).run()

            # save the corrected tree
            corrected_tree_path = os.path.join(self.saving_path, 'master_tree.corrected.nwk')
            phylo.to_file(self.tree, corrected_tree_path)
            self.tree_string = phylo.to_string(self.tree)

        self.taxaIndex, self.reverse = phylo.generate_taxa_index(self.tree, self.tax_filter, self.tax_mask)
        
        with open( self.saving_path + 'taxaIndex.pkl', 'wb') as taxout:
            taxout.write( pickle.dumps(self.taxaIndex))
        self.numperm = numperm
        ### if no weights are provided, generate them
        if treeweights is None:
            #generate aconfig_utilsll ones
            self.treeweights = hashutils.generate_treeweights(self.tree, self.taxaIndex, self.tax_filter, self.tax_mask)
        else:
            #load machine learning weights
            self.treeweights = treeweights
        tax_max = max(self.taxaIndex.values())+1
        wmg = WeightedMinHashGenerator(3*tax_max , sample_size = numperm , seed=1)
        with open( self.saving_path  + 'wmg.pkl', 'wb') as wmgout:
            wmgout.write( pickle.dumps(wmg))
        self.wmg = wmg

        print( '\nConfiguring pyham functions')

        print( 'taxfilter', self.tax_filter)
        print( 'taxmask', self.tax_mask)
        print( 'configuring pyham functions')
        print( 'swap ids', self.swap2taxcode)
        print( 'reformat names', self.reformat_names)
        print( 'use taxcodes', self.swap2taxcode)
        hamfunction = pyhamutils.get_ham_treemap_from_row
        hashfunction = hashutils.row2hash
        if slicesubhogs:
            hamfunction = pyhamutils.get_subhog_ham_treemaps_from_row
            hashfunction = hashutils.hash_trees_subhogs
        ### set up the pyHAM pipeline with different parameters depending on whether the input is an OMA hdf5 file or a list of orthoxml files
        print( 'lossonly', lossonly)
        print( 'duplonly', duplonly)

        self.dataset_nodes = None

        # remap taxfilter and taxmask 
        if taxfilter:
            self.tax_filter = [ self.idmapper[tax] for tax in taxfilter ]
            unacceptable_nodes = []
            for filterobj in self.tax_filter:
                filter_node = self.tree.search_nodes(name=filterobj)
                print('Found filter node:', filter_node)
                unacceptable_nodes.extend([node.name for node in filter_node[0].traverse()])
                
                #    print(f"Error searching for node with name: {filterobj}")
                #    print("Clade could not be excluded")
                #    continue
            # update dataset_nodes to exclude the filtered nodes
            self.dataset_nodes = [node.name for node in self.tree.traverse() if node.name not in unacceptable_nodes]

        if taxmask:
            #print(self.idmapper)
            self.tax_mask = self.idmapper[taxmask]
            ### get acceptable ids here:
            tax_mask_node = self.tree.search_nodes(name=self.tax_mask)
            if tax_mask_node:
                tax_mask_node = tax_mask_node[0]
                print(f"Found tax_mask_node: {tax_mask_node.name}")
                self.dataset_nodes = [node.name for node in tax_mask_node.traverse()]
            else:
                raise RuntimeError(f"No node found with name: {self.tax_mask}")
        
        print("Kept dataset nodes:", self.dataset_nodes)

        if self.h5OMA:
            self.HAM_PIPELINE = functools.partial( hamfunction, tree_string=self.tree_string ,  swap_ids=self.swap2taxcode , reformat_names = self.reformat_names , 
                                                  orthoXML_as_string = True, orthomapper = self.idmapper , levels = None,
                                                  limit_species = self.limit_species, limit_events = self.limit_events , dataset_nodes = self.dataset_nodes,
                                                  verbose = self.verbose) 
        else:
            self.HAM_PIPELINE = functools.partial( hamfunction, tree_string=self.tree_string ,  swap_ids=self.swap2taxcode  , 
                                                  orthoXML_as_string = False , reformat_names = self.reformat_names, 
                                                  orthomapper = self.idmapper , levels = None, limit_species = self.limit_species, limit_events = self.limit_events,
                                                  dataset_nodes = self.dataset_nodes, verbose = self.verbose)         
        ### set up the hash pipeline
        self.HASH_PIPELINE = functools.partial( hashfunction , taxaIndex=self.taxaIndex, treeweights=self.treeweights, wmg=wmg , lossonly = lossonly, duplonly = duplonly)
        print("\nSetting up input data reader")
        ### if the input is an OMA hdf5 file, set up the function to read the orthoxml data
        if self.h5OMA:
            self.READ_ORTHO = functools.partial(pyhamutils.get_orthoxml_oma, db_obj=self.db_obj)
            self.n_groups  = len(self.h5OMA.root.OrthoXML.Index)
            print( 'reading oma hdf5 with n groups:', self.n_groups)
        ### if the input is a list of orthoxml files, set up another function to read the orthoxml data
        elif self.fileglob:
            print('reading orthoxml files:' , len(self.fileglob))
            self.n_groups = len(self.fileglob)
        else:
            raise Exception( 'please specify an input file' )
        
        self.hashes_path = self.saving_path + 'hashes.h5'
        self.lshpath = self.saving_path + 'newlsh.pkl'
        self.lshforestpath = self.saving_path + 'newlshforest.pkl'
        self.mat_path = self.saving_path+ 'hogmat.h5'
        self.columns = len(self.taxaIndex)
        self.verbose = verbose
        print('done\n')

    def _get_species_names(self):
        """Read the DB or orthoxml to extract the list of species"""
        if self.h5OMA:
            field = "NCBITaxonId" if self.swap2taxcode else "SciName"
            values = self.h5OMA.root.Genome.read(field=field)
            return {
                value.decode("utf-8") if isinstance(value, bytes) else str(value)
                for value in values
            }
        else:
            names = set()
            for filename in self.fileglob:
                root = ET.parse(filename).getroot()
                for child in root:
                    if child.tag.rsplit("}", 1)[-1] == "species":
                        names.add(child.attrib["name"])
            return names

    def load_one(self, fam):
        #test function to try out the pipeline on one orthoxml
        ortho_fam = self.READ_ORTHO(fam)
        pyham_tree = self.HAM_PIPELINE([fam, ortho_fam])
        hog_matrix,weighted_hash = hashutils.hash_tree(pyham_tree , self.taxaIndex , self.treeweights , self.wmg)
        return ortho_fam , pyham_tree, weighted_hash,hog_matrix

    def generates_dataframes(self, size=100, minhog_size=10, maxhog_size=None):
        families = {}
        start = -1
        #hogsizetable = "/home/agavriil/Documents/venom_project/A_venom_analysis_tidy/2_profiling/1_oma_profiles/" \
        #"levels_oma_full_subhogs_250523/hogsize_table.csv"
        if self.h5OMA:
            self.groups  = self.h5OMA.root.OrthoXML.Index
            ### only Fam makes sense here, the rest are OMAmer related fields
            #print(self.h5OMA.root.OrthoXML.Index.colnames)
            self.rows = len(self.groups)
            for i, row in enumerate(tqdm.tqdm(self.groups)):
                if i > start:
                    #### family here is HOG ID minus the "HOG:E" prefix
                    fam = row[0]
                    ### testing only
                    #if fam != 794303 and fam!= 751313: #and fam != 712236 and fam != 708323:
                    #    continue
                    ortho_fam = self.READ_ORTHO(fam)
                    ### for older versions of OMA that do not have subhogIDs
                    orthoxmlbuilder = OrthoXMLBuilder()
                    ortho_fam = orthoxmlbuilder.add_loft_ids(ortho_fam)
                    hog_size = ortho_fam.count('<species name=')
                    #print(fam, hog_size)
                    #print(self.idmapper)
                    ### filtering for size already here (max HOG size and minimum species in HOG)
                    if (maxhog_size is None or hog_size < maxhog_size) and (minhog_size is None or hog_size > minhog_size):
                        families[fam] = {'ortho': ortho_fam}
                    if len(families) >= size:
                        pd_dataframe = pd.DataFrame.from_dict(families, orient='index')
                        pd_dataframe['Fam'] = pd_dataframe.index
                        yield pd_dataframe
                        families = {}
            if families:
                pd_dataframe = pd.DataFrame.from_dict(families, orient='index')
                pd_dataframe['Fam'] = pd_dataframe.index
                yield pd_dataframe

        elif self.fileglob:
            for i,file in enumerate(tqdm.tqdm(self.fileglob)):
                with open(file) as ortho:
                    #print("reading orthoxml file", file)
                    #oxml = ET.parse(ortho)
                    #ortho_fam = ET.tostring( next(oxml.iter()), encoding='utf8', method='xml' ).decode()
                    orthostr = ortho.read()
                    #print(os.path.basename(file),orthostr)
                hog_size = orthostr.count('<species name=')
                #print(os.path.basename(file),hog_size)
                if (maxhog_size is None or hog_size < maxhog_size) and (minhog_size is None or hog_size > minhog_size):
                    #print('fam',os.path.basename(file), 'saved')
                    families[i] = {'ortho': file}
                if len(families) >= size:
                    pd_dataframe = pd.DataFrame.from_dict(families, orient='index')
                    pd_dataframe['Fam'] = pd_dataframe.index
                    #print(pd_dataframe)
                    yield pd_dataframe
                    families = {}

            if families:
                pd_dataframe = pd.DataFrame.from_dict(families, orient='index')
                pd_dataframe['Fam'] = pd_dataframe.index
                yield pd_dataframe


    #def worker(self, i, work_queue, result_queue):
    #    logger.debug('Starting worker #%d', i)

    def worker(self, i, q, retq):
        try:
            if self.verbose == True:
                print('worker init ' + str(i))
            while True:
                ### get dataframe from queue (generates_dataframes)
                df = q.get()
                #print(df.head()) # gives error if df is None
                if df is not None :
                    if self.slicesubhogs is False:
                        df['tree'] = df[['Fam', 'ortho']].apply(self.HAM_PIPELINE, axis=1)
                        #add a dictionary of results with subhogs { fam_sub1: { 'tree':tp , 'Fam':fam }  , fam_sub2: { 'tree':tp , 'Fam':fam } , ... }
                        #returned_df = pd.DataFrame.from_dict(df['tree'].to_dict(), orient='index')
                        #merge with pandas on right e.g. df.merge( returned_df , on = 'Fam' , how = 'right' )
                        #print(df.head())
                        df[['hash','rows']] = df[['Fam', 'tree']].apply(self.HASH_PIPELINE, axis=1)
                        ### check for empty dfs before putting them in the queue
                        if df[['Fam', 'hash']].empty:
                            return
                        ### put data in queue to be saved
                        if self.fileglob:
                            retq.put(df[['Fam', 'hash', 'ortho']])
                        else:
                            retq.put(df[['Fam', 'hash']])
                    ### slice subhogs case
                    else:
                        #print(df)
                        #print(df.size)
                        ### original one liner
                        #df['tree_dicts'] = df[['Fam', 'ortho']].apply(self.HAM_PIPELINE, axis=1)
                        ### for debugging
                        # Iterate over each row in the DataFrame
                        tree_dicts = []
                        #print(f"Memory before HAM: {process.memory_info().rss / 1024 / 1024 / 1024:.2f} GB")
                        for index, row in df[['Fam', 'ortho']].iterrows():
                            try:
                                # Apply the HAM_PIPELINE function to the current row
                                result = self.HAM_PIPELINE(row)
                                tree_dicts.append(result)
                                #print(f"Debug: Successfully processed row {index} with Fam={row['Fam']}")
                            except Exception as e:
                                # Handle and log any errors
                                print(f"Error processing row {index} with Fam={row['Fam']}: {e}")
                                # Append None for rows that failed - necessary or raises alueError: Length of values (65) does not match length of index (101)
                                tree_dicts.append(None)  
                        #print(f"Memory after HAM: {process.memory_info().rss / 1024 / 1024 / 1024:.2f} GB")
                        #df['tree_dicts'] = tree_dicts
                        #print(tree_dicts)
                        # Assign the results back to the DataFrame
                        df['tree_dicts'] = tree_dicts
                        df['hash_dicts'] = df[['Fam', 'tree_dicts']].apply(self.HASH_PIPELINE, axis=1)
                        #print(f"Memory after HASH: {process.memory_info().rss / 1024 / 1024 / 1024:.2f} GB")
                        #print(df)
                        #print(df.tree_dicts.iloc[0])
                        # Filter out rows with empty hash_dicts to save time
                        df = df[df['hash_dicts'].apply(bool)]
                        ### check again if df is empty:
                        if df.empty:
                            return
                        #print(df.tree_dicts.iloc[0])
                        newdf ={}
                        #print("worker df columns",df.columns)
                        #for col in df.columns:
                        #    print(df[[col]].head())
                        for i,row in df.iterrows():
                            for subhog in row['hash_dicts']:
                                if self.fileglob:
                                    newdf[ ( row['Fam'] ,  subhog ) ] = { 'tree': row['tree_dicts'][subhog] , 'hash': row['hash_dicts'][subhog][1] 
                                                                    , 'ortho': row['ortho'] }
                                ### don't save orthoxml strings if OMA
                                else:
                                    newdf[ ( row['Fam'] ,  subhog ) ] = { 'tree': row['tree_dicts'][subhog] , 'hash': row['hash_dicts'][subhog][1],
                                    'ortho':''}
                        newdf = pd.DataFrame.from_dict(newdf, orient='index')
                        ### this here is what prints empty for OMA run (seemingly cause of the Toxicofera mask)!!!!!!!!!!!!
                        #print(newdf)
                        #print(f"Memory before retq: {process.memory_info().rss / 1024 / 1024 / 1024:.2f} GB")
                        retq.put(newdf)
                else:
                    if self.verbose == True:
                        print('Worker done' + str(i))
                    break
        except Exception as e:
            import traceback
            print('Worker error', file=sys.stderr)
            print(f"Error in worker process: {traceback.format_exc()}", file=sys.stderr)
            
    
    def _make_tax_str(self):
        taxstr = ""
        if self.tax_filter is None:
            taxstr = "NoFilter"
        if self.tax_mask is None:
            taxstr += "NoMask"
        else:
            taxstr = str(self.tax_filter)
        return taxstr

    def _index_and_save(self, forest):
        logger.debug("Saving forest at: %.2f", t.time() - self.start_time)

        forest.index()
        with open(self.lshforestpath, "wb") as forest_out:
            forest_out.write(pickle.dumps(forest, -1))

        logger.debug("Save done at: %.2f", t.time() - self.start_time)

    def _save_family_mapping(self, frames):
        if frames:
            mapping = pd.concat(frames)
        else:
            mapping = pd.DataFrame(columns=['Fam', 'ortho'])
        mapping.to_csv(self.saving_path + 'fam2orthoxml.csv')


    def saver(self, result_queue):
        save_start = t.time()
        global_time = t.time()
        self.start_time = global_time

        chunk_size = 100
        count = 0
        last_reported_count = 0
        forest = MinHashLSHForest(num_perm=self.numperm)
        mapping_frames = []
        taxstr = self._make_tax_str()

        with h5py.File(self.hashes_path, 'w', libver='latest') as h5hashes:
            hash_width = 2 * self.numperm
            dataset = h5hashes.create_dataset(
                taxstr,
                (chunk_size, hash_width),
                maxshape=(None, hash_width),
                dtype="int32",
            )
            h5hashes.flush()
            logger.debug('Creating dataset filtered at taxonomic level: %s', taxstr)

            done = False
            while not done:
                this_dataframe = result_queue.get()
                if this_dataframe is not None:
                    if not this_dataframe.empty:
                        hashes = this_dataframe['hash'].to_dict()
                        #print(str(this_dataframe.Fam.max())+ 'fam num')
                        #print(str(count) + ' done')

                        hashes = { fam:hashes[fam] for fam in hashes if hashes[fam]}
                        for fam in hashes:
                            forest.add(str(fam), hashes[fam])

                        for fam in hashes:
                            # if all rows are filled, allocate the next chunk
                            if dataset.shape[0] <= fam:
                                dataset.resize(fam + chunk_size, axis=0)

                            dataset[fam, :] = hashes[fam].hashvalues.ravel()
                            count += 1

                        if count - last_reported_count >= 1000:
                            # maybe we don't need this as this disrupts tqdm output
                            #logger.info('Saver has written %d hashes', count)
                            last_reported_count = count

                        if self.fileglob:
                            mapping_frames.append(this_dataframe[['Fam', 'ortho']])

                        if hashes and t.time() - save_start > 200:
                            h5hashes.flush()
                            self._index_and_save(forest)
                            logger.debug('Testing forest')
                            logger.debug(forest.query(hashes[fam], k=10))

                            if self.fileglob:
                                self._save_family_mapping(mapping_frames)
                            save_start = t.time()
                    else:
                        print(this_dataframe)
                else:

                    logger.info('Wrapping up the run')

                    self._index_and_save(forest)

                    h5hashes.flush()
                    if self.fileglob:
                        self._save_family_mapping(mapping_frames)

                    logger.info('Saver wrote %d hashes', count)
                    done = True


    @staticmethod
    def _raise_if_failed(processes, required_alive=()):
        """Raise if a child process failed"""
        failed = [
            process for process in processes
            if process.exitcode not in (None, 0)
        ]
        if failed:
            details = ", ".join(
                f"{process.name} (exit code {process.exitcode})"
                for process in failed
            )
            raise RuntimeError(f"Child process failure: {details}")

        stopped = [
            process for process in required_alive
            if process.exitcode is not None
        ]
        if stopped:
            names = ", ".join(process.name for process in stopped)
            raise RuntimeError(f"Child process exited unexpectedly: {names}")

    @classmethod
    def _safe_put(cls, destination, value, processes, required_alive=()):
        """
        Enqueue data value with process/queue checks:
            - if any of the workers failed already, don't send
            - if queue is full, try again
        """
        while True:
            cls._raise_if_failed(processes, required_alive)

            try:
                destination.put(value, timeout=0.2)
                return
            except queue.Full:
                pass

    @classmethod
    def _wait_for_all(cls, processes, monitored_processes=()):
        """
        Joins a list of processes; raises if any of them
        failed or if the monitored processes failed
        """
        alive = list(processes)
        all_processes = [*alive, *monitored_processes]

        while alive:
            for process in alive:
                process.join(timeout=0)

            cls._raise_if_failed(all_processes)

            alive = [process for process in alive if process.is_alive()]
            if alive:
                time.sleep(0.1)

    @staticmethod
    def _terminate_all(processes):
        for process in processes:
            if process.is_alive():
                process.terminate()

        for process in processes:
            process.join()


    def _run_parallel(self, threads, data_generator):
        if threads < 1:
            raise ValueError("threads must be at least 1")

        queue_size = max(2, threads * 2)
        work_queue = mp.Queue(maxsize=queue_size)
        result_queue = mp.Queue(maxsize=queue_size)
        saver = mp.Process(name="saver", target=self.saver, args=(result_queue,))
        workers = [
            mp.Process(name=f"worker-{i}", target=self.worker, args=(i, work_queue, result_queue))
            for i in range(threads)
        ]
        processes = [saver, *workers]
        started_processes = []
        completed = False

        try:
            # start saver
            saver.start()
            started_processes.append(saver)
            logger.info('Started one saver process')

            # start workers
            for process in workers:
                process.start()
                started_processes.append(process)
            logger.info('Started %d workers', threads)

            # send data batches to workers
            for data in data_generator:
                self._safe_put(work_queue, data, processes, required_alive=processes)

            # send sentinels (end of work signals) to all workers
            for _ in workers:
                self._safe_put(work_queue,None, processes, required_alive=(saver,))
            logger.info('Input queued; waiting for workers to finish')

            # the main join -- wait for workers to be done
            self._wait_for_all(workers, monitored_processes=(saver,))

            logger.info('Workers finished; finalizing saver')
            # send the sentinel to the saver
            self._safe_put(result_queue,None,(saver,), required_alive=(saver,))

            # wait for the saver
            self._wait_for_all((saver,))
            completed = True

            if self.slicesubhogs:
                csv_path = os.path.join(self.saving_path, 'fam2orthoxml.csv')
                if os.path.exists(csv_path):
                    df = pd.read_csv(csv_path, index_col=[0, 1])  # Read first two columns as index
                    #print(df.head())
                    #print(df.size)
                    df.index.set_names(['fam', 'subhog_id'], inplace=True)  # Set correct names
                    df.to_csv(csv_path)  # Overwrite with corrected index names
                    #print(df.size)

        finally:
            # if error
            if not completed:
                self._terminate_all(started_processes)

            # finalize queues
            for process_queue in (work_queue, result_queue):
                if not completed:
                    # if failed, the join below is ignored
                    process_queue.cancel_join_thread()

                process_queue.close()
                if completed:
                    process_queue.join_thread()

        logger.info('Pipeline complete')

    def safe_ham(self, row):
            try:
                return self.HAM_PIPELINE(row)
            except Exception as e:
                if self.verbose:
                    print(f"Error in HAM_PIPELINE for Fam={row['Fam']}: {e}")
                return None

    def worker_single(self, i, data, retq, matq, l):
        if self.verbose:
            print('worker init ' + str(i))

        if data is not None:
            df = data
            #print(df.head())
            if self.slicesubhogs is False:
                df['tree'] = df[['Fam', 'ortho']].apply(self.safe_ham, axis=1)
                # Drop rows where tree is None (error in HAM_PIPELINE)
                df = df[df['tree'].notnull()]
                #print(f'Generated tree in worker_single: {df["tree"]}')
                #print("df fam tree\n", df[['Fam', 'tree']])  # Debugging
                df[['hash', 'rows']] = df[['Fam', 'tree']].apply(self.HASH_PIPELINE, axis=1)
                #print(f'Generated hash and rows in worker_single: {df[["hash", "rows"]]}')
                if self.fileglob:
                    return df[['Fam', 'hash', 'ortho']]
                else:
                    return df[['Fam', 'hash']]
            else:

                df['tree_dicts'] = df[['Fam', 'ortho']].apply(self.safe_ham, axis=1)
                # Drop rows where tree is None (error in HAM_PIPELINE)
                df = df[df['tree_dicts'].notnull()]
                #print(df)
                df['hash_dicts'] = df[['Fam', 'tree_dicts']].apply(self.HASH_PIPELINE, axis=1)
                
                newdf ={}
                for i,row in df.iterrows():
                    for subhog in row['hash_dicts']:
                        if self.fileglob:
                            newdf[ ( row['Fam'] ,  subhog ) ] = { 'tree': row['tree_dicts'][subhog] , 'hash': row['hash_dicts'][subhog][1] 
                                                             , 'ortho': row['ortho'] }
                        else:
                            newdf[ ( row['Fam'] ,  subhog ) ] = { 'tree': row['tree_dicts'][subhog] , 'hash': row['hash_dicts'][subhog][1] }
                newdf = pd.DataFrame.from_dict(newdf, orient='index')
                #print(newdf)
                empty_rows = df[df['tree_dicts'].apply(lambda x: len(x) == 0)]
                if self.verbose:
                    print(f'{empty_rows.shape[0]} famililes failed the filter check')
                    print(empty_rows)

                return newdf



    ### SAVER LEVELS
    #
    ### saving fams to orthoxml mapping to a csv file (fam2orthoxml.csv)
    ### creating error file to log any errors (errors.txt)
    ### saving the MinHashLSHForest to a pickle file (newlshforest.pkl)
    ### saving the MinHashes to an hdf5 file (hashes.h5)
    # def saver(self, retq):
    #     try:
    #         print_start = t.time()
    #         save_start = t.time()
    #         global_time = t.time()
    #         chunk_size = 100
    #         count = 0
    #         forest = MinHashLSHForest(num_perm=self.numperm)
    #         taxstr = ''
    #         savedf = None
    #         total_subfam_ids = []
    #         totals_subfams = 0
    #         if self.tax_filter is None:
    #             taxstr = 'NoFilter'
    #         if self.tax_mask is None:
    #             taxstr+= 'NoMask'
    #         else:
    #             taxstr = str(self.tax_filter)
    #         self.errorfile = self.saving_path + 'errors.txt'
    #         with open(self.errorfile, 'w') as hashes_error_files:
    #             with h5py.File(self.hashes_path, 'w', libver='latest') as h5hashes:
                
    #                 hash_width = 2 * self.numperm
    #                 dataset = h5hashes.create_dataset(
    #                     taxstr,
    #                     (chunk_size, hash_width),
    #                     maxshape=(None, hash_width),
    #                     dtype="int32",
    #                 )
    #                 h5hashes.flush()
    #                 logger.debug('Creating dataset filtered at taxonomic level: %s', taxstr)

    #                 while True:
    #                     #print("Debug: Checking retq queue before processing:")
    #                     #print(retq.qsize())
    #                     this_dataframe = retq.get()
    #                     if this_dataframe is not None:
    #                         if not this_dataframe.empty:
    #                             hashes = this_dataframe['hash'].to_dict()
    #                             #print(str(this_dataframe.Fam.max())+ 'fam num')
    #                             #print(str(count) + ' done')
    #                             ### remove empty hashes
    #                             hashes = {fam:hashes[fam]  for fam in hashes if hashes[fam] }
    #                             if self.verbose == True:
    #                                 print(f'Non empty hashes: {len(hashes)}')
    #                                 print(this_dataframe)
    #                             ### handle slicesubhogs
    #                             if self.slicesubhogs:
    #                                 subfam_ids_list = this_dataframe.index.to_list()
    #                                 total_subfam_ids.extend(subfam_ids_list)
    #                                 nsubfams = len(hashes)
    #                                 # Resize if necessary
    #                                 if h5hashes[taxstr].shape[0] < totals_subfams + nsubfams + chunk_size:
    #                                     example = list(hashes.keys())[0]
    #                                     h5hashes[taxstr].resize((totals_subfams + nsubfams + chunk_size, len(hashes[example].hashvalues.ravel())))
    #                                 # Store hashes
    #                                 h5hashes[taxstr][totals_subfams:totals_subfams + nsubfams, :] = [hashes[fam].hashvalues.ravel() for fam in hashes]

    #                                 # Update forest
    #                                 [forest.add(str(fam[0]) + '_' + str(fam[1]), hashes[fam]) for fam in hashes]
                                
    #                                 totals_subfams += nsubfams

    #                                 #if self.fileglob or self.slicesubhogs:
    #                                 if savedf is None:
    #                                     #df_cols = this_dataframe.columns
    #                                     savedf = this_dataframe[['ortho']]
    #                                 else:
    #                                     savedf = pd.concat([savedf, this_dataframe[['ortho']]])

    #                             ### standard processing
    #                             else:
    #                                 [forest.add(str(fam), hashes[fam]) for fam in hashes]
    #                                 for fam in hashes:
    #                                     if len(h5hashes[taxstr]) < fam + 10:
    #                                         h5hashes[taxstr].resize((fam + chunk_size, len(hashes[fam].hashvalues.ravel())))
    #                                     h5hashes[taxstr][fam, :] = hashes[fam].hashvalues.ravel()
    #                                 if self.fileglob:# or self.slicesubhogs:
    #                                     # Addition to ensure savedf is properly updated - 13.05.25
    #                                     if not this_dataframe.empty:
    #                                         if self.verbose:
    #                                             print(this_dataframe)
    #                                         if savedf is None:
    #                                             savedf = this_dataframe[['Fam', 'ortho']]
    #                                         else:
    #                                             savedf = pd.concat([savedf, this_dataframe[['Fam', 'ortho']]])
    #                                         if self.verbose:
    #                                             print(savedf)
                                    
    #                             # Save every 200 seconds
    #                             if t.time() - save_start > 200:
    #                                 print('Saving at:', t.time() - global_time)
    #                                 forest.index()
    #                                 print( 'testing forest' )
    #                                 testfam = list(hashes.keys())[0]
    #                                 print(forest.query( hashes[testfam] , k = 10 ) )
    #                                 h5flush()
    #                                 with open(self.lshforestpath, 'wb') as forestout:
    #                                     forestout.write(pickle.dumps(forest, -1))
    #                                 if self.fileglob or self.slicesubhogs:
    #                                     #print('Saving fam-to-orthoxml mapping')
    #                                     savedf.to_csv(os.path.join(self.saving_path, 'fam2orthoxml.csv'))
    #                                 save_start = t.time()

    #                     # wrap up
    #                     else:
    #                         print('\nwrapping up the run')
    #                         print('saving at :' , t.time() - global_time )
    #                         forest.index()
    #                         with open(self.lshforestpath , 'wb') as forestout:
    #                             forestout.write(pickle.dumps(forest, -1))
    #                         h5flush()
    #                         ### make sure fam2orthoxml is saved no matter what
    #                         if self.slicesubhogs:
    #                             if savedf is not None and not savedf.empty:
    #                                 print('saving orthoxml to fam mapping')
    #                                 print(savedf.head())
    #                                 savedf.to_csv(os.path.join(self.saving_path, 'fam2orthoxml.csv'))
    #                             elif self.h5OMA and savedf is not None:
    #                                 # Save a mapping with just Fam and subhog_id (index of savedf)
    #                                 mapping = pd.DataFrame(savedf.index.tolist(), columns=['Fam', 'subhog_id'])
    #                                 mapping.to_csv(os.path.join(self.saving_path, 'fam2orthoxml.csv'), index=False)
    #                             else:
    #                                 print("savedf is empty?")
    #                                 print('Warning: No fam-to-orthoxml mapping found.')
    #                         print('DONE SAVER')
    #                         break
    #     except Exception as e:
    #         import traceback
    #         print('Worker error')
    #         print(f"Error in worker process: {traceback.format_exc()}")


    ### multithreaded pipeline to compile the LSH forest database of MinHashes
    #def run_pipeline(self , threads):
        # print( 'run w n threads:', threads)
        # functype_dict = {'worker': (self.worker, threads , True), 'updater': (self.saver, 1, False),
        #                  'matrix_updater': (self.matrix_updater, 0, False) }
        # ### function to manage parallel processing
        # def mp_with_timeout(functypes, data_generator):
        #     # variables to store processes
        #     work_processes = {} 
        #     update_processes = {}
        #     # lock to synchronize access to shared resources
        #     lock = mp.Lock()
        #     cores = mp.cpu_count()
        #     # objects to send data between processes
        #     q = mp.Queue(maxsize=cores * 10)
        #     retq = mp.Queue(maxsize=cores * 10)
        #     matq = mp.Queue(maxsize=cores * 10)
        #     work_processes = {}
        #     error_queue = mp.Queue() 
        #     print('start workers')
        #     for key in functypes:
        #         worker_function, number_workers, joinval = functypes[key]
        #         work_processes[key] = []
        #         for i in range(int(number_workers)):
        #             t = mp.Process(target=worker_function, args=(i, q, retq, matq, lock ))
        #             t.daemon = True
        #             work_processes[key].append(t)
        #     # starting the processes for each worker in parallel
        #     for key in work_processes:
        #         for process in work_processes[key]:
        #             process.start()
        #     # feeding data to workers - data from generator are pushed into the queue q for the workers
        #     for data in tqdm.tqdm(data_generator):
        #         q.put(data)
            
        #     print('done spooling data')
        #     # stopping workers - after all data are processed, None is pushed into the queue to 
        #     # signal the workers to stop
        #     for key in work_processes:
        #         for i in range(2):
        #             for _ in work_processes[key]:
        #                 q.put(None)
        #     print('joining processes')
        #     # if joinval is True, the main process waits for the workers to finish
        #     for key in work_processes:
        #         worker_function, number_workers , joinval = functypes[key]
        #         if joinval == True:
        #             for process in work_processes[key]:
        #                 process.join()
        #     # processing data for non-joinable workers (updater and matrix updater) - pushing their
        #     # results into retq and matq for processing
        #     for key in work_processes:
        #         worker_function, number_workers, joinval = functypes[key]
        #         if joinval == False:
        #             for _ in work_processes[key]:
        #                 retq.put(None)
        #                 matq.put(None)
        #     # final joining - waiting for all non-joinable processes to finish before proceeding
        #     for key in work_processes:
        #         worker_function, number_workers , joinval = functypes[key]
        #         if joinval == False:
        #             for process in work_processes[key]:
        #                 process.join()
        #     # Check error queue
        #     while not error_queue.empty():
        #         print("ERROR in worker process:")
        #         print(error_queue.get())
        #     # garbage collection
        #     gc.collect()
        #     print('DONE!')
        # ### run the pipeline using the functions defined above
        # limit_species = self.limit_species
        # mp_with_timeout(functypes=functype_dict, data_generator=self.generates_dataframes(size=100, minhog_size=limit_species)) #original was 100
        # ### fix labels for slicesubhogs
        # #'''
        # if self.slicesubhogs:
        #     csv_path = os.path.join(self.saving_path, 'fam2orthoxml.csv')
        #     if os.path.exists(csv_path):
        #         df = pd.read_csv(csv_path, index_col=[0, 1])  # Read first two columns as index
        #         #print(df.head())
        #         #print(df.size)
        #         df.index.set_names(['fam', 'subhog_id'], inplace=True)  # Set correct names
        #         df.to_csv(csv_path)  # Overwrite with corrected index names
        #         #print(df.size)
        # #'''
        # ### return the paths to the output files
        # return self.hashes_path, self.lshforestpath , self.mat_path
    
    def run_pipeline(self , threads):
        logger.info('Running with %d threads', threads)
        data_generator = self.generates_dataframes(size=100, minhog_size=self.limit_species)
        self._run_parallel(threads=threads, data_generator=data_generator)
        return self.hashes_path, self.lshforestpath, self.mat_path
    

    def run_pipeline_single(self):
        print('\nRun single-threaded pipeline')
        limit_species = self.limit_species
        # Use the correct data generator
        data_generator = self.generates_dataframes(size=100, minhog_size=limit_species)
        
        # Initialize minimal required variables
        chunk_size = 100
        forest = MinHashLSHForest(num_perm=self.numperm)
        savedf = None
        total_subfam_ids = []
        totals_subfams = 0
        global_time = t.time()
        
        # Set taxstr consistently with saver function
        if self.tax_filter is None:
            taxstr = 'NoFilter'
        if self.tax_mask is None:
            taxstr += 'NoMask'
        else:
            taxstr = str(self.tax_filter)
            
        # Open files for writing, matching saver's file handling
        with open(self.errorfile, 'w') as hashes_error_files:
            with h5py.File(self.hashes_path, 'w', libver='latest') as h5hashes:
                # Initialize h5 dataset
                if taxstr not in h5hashes.keys():
                    if self.verbose:
                        print('creating dataset')
                        print('filtered at taxonomic level: '+taxstr)
                    h5hashes.create_dataset(taxstr, (chunk_size, 0), maxshape=(None, None), dtype='int32')
                
                this_dataframe = None
                for i, data in enumerate(tqdm.tqdm(data_generator)):
                    if data is not None:
                        retdf = self.worker_single(i, data, [], [], None)  # Empty lists for retq/matq since not using queues
                        
                        if self.slicesubhogs is False:
                            # Handle non-sliced case
                            if self.fileglob:
                                if savedf is None:
                                    savedf = retdf[['Fam', 'ortho']]
                                else:
                                    savedf = pd.concat([savedf, retdf[['Fam', 'ortho']]])
                            
                            if this_dataframe is None:
                                this_dataframe = retdf
                            elif this_dataframe is not None:
                                this_dataframe = pd.concat([this_dataframe, retdf])
                            
                            if len(this_dataframe) > 10:
                                hashes = this_dataframe['hash'].to_dict()
                                hashes = {fam:hashes[fam] for fam in hashes if hashes[fam]}
                                
                                # Add to forest and resize h5 if needed
                                for fam in hashes:
                                    if len(h5hashes[taxstr]) < fam + chunk_size:
                                        h5hashes[taxstr].resize((fam + chunk_size, len(hashes[fam].hashvalues.ravel())))
                                    forest.add(str(fam), hashes[fam])
                                    h5hashes[taxstr][fam, :] = hashes[fam].hashvalues.ravel()
                                
                                this_dataframe = None
                                
                        else:
                            # Handle sliced case 
                            if self.fileglob:
                                if savedf is None:
                                    savedf = retdf[['ortho']]
                                else:
                                    savedf = pd.concat([savedf, retdf[['ortho']]])
                            
                            if this_dataframe is None:
                                this_dataframe = retdf
                            elif this_dataframe is not None:
                                this_dataframe = pd.concat([this_dataframe, retdf])
                            
                            if len(this_dataframe) > 10:
                                hashes = this_dataframe['hash'].to_dict()
                                nsubfams = len(hashes)
                                
                                # Resize if needed
                                if h5hashes[taxstr].shape[0] < totals_subfams + nsubfams + chunk_size:
                                    example = list(hashes.keys())[0]
                                    h5hashes[taxstr].resize((totals_subfams + nsubfams + chunk_size, 
                                                           len(hashes[example].hashvalues.ravel())))
                                
                                subfam_ids_list = retdf.index.to_list()
                                total_subfam_ids.extend(subfam_ids_list)
                                
                                # Store hashes
                                h5hashes[taxstr][totals_subfams:totals_subfams + nsubfams, :] = \
                                    [hashes[fam].hashvalues.ravel() for fam in subfam_ids_list]
                                
                                totals_subfams += nsubfams
                                
                                # Add to forest
                                hashes = {fam:hashes[fam] for fam in hashes if hashes[fam]}
                                for fam in hashes:
                                    forest.add(str(fam[0]) + '_' + str(fam[1]), hashes[fam])
                                
                                this_dataframe = None

                # Handle final dataframe if exists
                if this_dataframe is not None:
                    if self.slicesubhogs:
                        # Handle remaining sliced data
                        hashes = this_dataframe['hash'].to_dict() 
                        nsubfams = len(hashes)
                        if h5hashes[taxstr].shape[0] < totals_subfams + nsubfams + chunk_size:
                            example = list(hashes.keys())[0]
                            h5hashes[taxstr].resize((totals_subfams + nsubfams + chunk_size,
                                                   len(hashes[example].hashvalues.ravel())))
                        
                        h5hashes[taxstr][totals_subfams:totals_subfams + nsubfams, :] = \
                            [hashes[fam].hashvalues.ravel() for fam in hashes]
                        totals_subfams += nsubfams
                        
                        hashes = {fam:hashes[fam] for fam in hashes if hashes[fam]}
                        for fam in hashes:
                            forest.add(str(fam[0]) + '_' + str(fam[1]), hashes[fam])
                    
                    else:
                        # Handle remaining non-sliced data
                        hashes = this_dataframe['hash'].to_dict()
                        hashes = {fam:hashes[fam] for fam in hashes if hashes[fam]}
                        for fam in hashes:
                            if len(h5hashes[taxstr]) < fam + chunk_size:
                                h5hashes[taxstr].resize((fam + chunk_size, len(hashes[fam].hashvalues.ravel())))
                            forest.add(str(fam), hashes[fam])
                            h5hashes[taxstr][fam, :] = hashes[fam].hashvalues.ravel()

                # Final processing
                print('wrapping up the run')
                print('saving at:', t.time() - global_time)
                
                # Index and save forest
                forest.index()
                with open(self.lshforestpath, 'wb') as forestout:
                    forestout.write(pickle.dumps(forest, -1))
                
                # Save mappings if needed
                if savedf is not None:
                    if self.slicesubhogs:
                        savedf.index.set_names(['fam','subhog_id'], inplace=True)
                        savedf.reset_index(inplace=True)
                        savedf.index = range(len(savedf))
                    
                    if self.fileglob:
                        print('saving orthoxml to fam mapping')
                        savedf.to_csv(os.path.join(self.saving_path, 'fam2orthoxml.csv'))

        print('done single-threaded pipeline')

                

    

def main():
    parser = argparse.ArgumentParser()
    parser = argparse.ArgumentParser(prog="hogprof")
    parser.add_argument('--version', action='version',
                        version=f'%(prog)s {__version__}')
    parser.add_argument('--taxweights', help='load optimised weights from keras model',type = str)
    parser.add_argument('--taxmask', help='consider only one branch (e.g. Sauria)',type = str)
    parser.add_argument('--taxfilter', help='remove these taxa' , type = str)
    parser.add_argument('--outpath', help='name of the db (output folder where all files will be created)', type = str)
    parser.add_argument('--dbtype', help='preconfigured taxonomic ranges' , type = str)
    parser.add_argument('--OMA', help='use oma data ' , type = str)
    parser.add_argument('--OrthoGlob', help='a glob expression for orthoxml files ' , type = str)
    parser.add_argument('--tarfile', help='use tarfile with orthoxml data ' , type = str)
    parser.add_argument('--nperm', help='number of hash functions to use when constructing profiles' , type = int)
    parser.add_argument('--mastertree', help='master taxonomic tree. nodes should correspond to orthoxml' , type = str)
    
    parser.add_argument('--nthreads', help='nthreads for multiprocessing' , type = int)
    parser.add_argument('--lossonly', help='only compile loss events' , type = bool)
    parser.add_argument('--duplonly', help='only compile duplication events' , type = bool)
    parser.add_argument('--taxcodes', help='use taxid info in HOGs' , type = bool)
    parser.add_argument('--verbose', help='print verbose output', action='store_true')
    parser.add_argument('--reformat_names', help='try to correct broken species trees by replacing all names with numbers.', action='store_true')
    parser.add_argument('--slicesubhogs', help='slice subhogs' , type = bool, default=False)
    parser.add_argument('--specieslim', help='minimum number of species in a subhog' , type = int, default=10)
    parser.add_argument('--eventslim', help='minimum number of events (loss/duplication) in a subhog' , type = int, default=0)
    
    dbdict = {
        'all': { 'taxfilter': None , 'taxmask': None },
        'plants': { 'taxfilter': None , 'taxmask': 33090 },
        'archaea':{ 'taxfilter': None , 'taxmask': 2157 },
        'bacteria':{ 'taxfilter': None , 'taxmask': 2 },
        'eukarya':{ 'taxfilter': None , 'taxmask': 2759 },
        'protists':{ 'taxfilter': [2 , 2157 , 33090 , 4751, 33208] , 'taxmask':None },
        'fungi':{ 'taxfilter': None , 'taxmask': 4751 },
        'metazoa':{ 'taxfilter': None , 'taxmask': 33208 },
        'vertebrates':{ 'taxfilter': None , 'taxmask': 7742 },
    }

    taxfilter = None
    taxmask = None
    omafile = None

    args = vars(parser.parse_args(sys.argv[1:]))
    print("\nReading arguments")

    orthoglob = None
    if 'OrthoGlob' in args:
        if args['OrthoGlob']:
            print("Using orthoxml files from:", args['OrthoGlob'])
            orthoglob = glob.glob(args['OrthoGlob'])

    if 'outpath' in args:
        dbname = args['outpath']
        if not dbname.endswith('/'):
            dbname += '/'
    else:
        raise Exception(' please give your profile an output path with the --outpath argument ')
    if args['dbtype']:
        taxfilter = dbdict[args['dbtype']]['taxfilter']
        taxmask = dbdict[args['dbtype']]['taxmask']
    if args['taxmask']:
        taxmask = args['taxmask']
    
    if args['taxfilter']:
        taxfilter = args['taxfilter']

    if args['nperm']:
        nperm = int(args['nperm'])
    else:
        nperm = 256
    if args['OMA']:
        omafile = args['OMA']
    elif args['tarfile']:
        omafile = args['tarfile']
    elif orthoglob is not None:
        fileglob = orthoglob
    else:
        raise Exception(' please specify input data ')
    if args['lossonly']:
        lossonly = args['lossonly']
    else:
        lossonly = False
    if args['duplonly']:
        duplonly = args['duplonly']
    else:
        duplonly = False

    if args['taxcodes']==True:
        taxcodes = True
    else:
        taxcodes = False

    _args = parser.parse_args()
    verbose = _args.verbose

    if args['reformat_names']:
        reformat_names = True
    else:
        reformat_names = False

    threads = 4
    if args['nthreads']:
        threads = args['nthreads']
    if args['taxweights']:
        from keras.models import model_from_json
        json_file = open(  args['taxweights']+ '.json', 'r')
        loaded_model_json = json_file.read()
        json_file.close()
        model = model_from_json(loaded_model_json)
        # load weights into new model
        model.load_weights(  args['taxweights']+".h5")
        print("Loaded model from disk")
        weights = model.get_weights()[0]
        weights += 10 ** -10
    else:
        weights = None
    if args['mastertree']:
        mastertree = Path(args['mastertree'])
    else:
        mastertree=None

    # set DEBUG log level if --verbose
    log_level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(level=logging.INFO)

    # Pyham spams INFO-level messages like crazy. Suppress
    logging.getLogger("pyham").setLevel(logging.WARNING)
    logging.getLogger().setLevel(log_level)

    start = time.time()
    if omafile:
        with open_file( omafile , mode="r") as h5_oma:
            lsh_builder = LSHBuilder(h5_oma = h5_oma,  fileglob=orthoglob ,saving_name=dbname , numperm = nperm ,
            treeweights= weights , taxfilter = taxfilter, taxmask=taxmask , masterTree =mastertree , 
            lossonly = lossonly , duplonly = duplonly , use_taxcodes = taxcodes , reformat_names=reformat_names, 
            verbose=verbose, slicesubhogs=args['slicesubhogs'], limit_species=args['specieslim'], limit_events=args['eventslim'])
            lsh_builder.run_pipeline(threads)
            #lsh_builder.run_pipeline_single() # made for local tests. ignore

    else:
        lsh_builder = LSHBuilder(h5_oma = None,  fileglob=orthoglob ,saving_name=dbname , numperm = nperm ,
        treeweights= weights , taxfilter = taxfilter, taxmask=taxmask ,
          masterTree =mastertree , lossonly = lossonly , duplonly = duplonly , use_taxcodes = taxcodes , 
          reformat_names=reformat_names, verbose=verbose, slicesubhogs=args['slicesubhogs'], limit_species=args['specieslim'], 
          limit_events=args['eventslim'])
        lsh_builder.run_pipeline(threads)

    logger.info("Done in %.2fs", time.time() - start)


if __name__ == '__main__':
    main()
