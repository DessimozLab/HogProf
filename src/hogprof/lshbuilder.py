import functools
import glob
import logging
import multiprocessing as mp
import os
import pickle
import queue
import random
import time
import time as t
from datetime import datetime
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
from datasketch import MinHashLSHForest, WeightedMinHashGenerator
from pyoma.browser import db
from tables import open_file

from hogprof.cli import DB_PRESETS, track_progress
from hogprof.utils import hashutils, phylo, pyhamutils
from hogprof.utils import orthoxml as oxml

logger = logging.getLogger(__name__)

random.seed(0)
np.random.seed(0)


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

    def __init__(self, 
                 h5_oma=None, fileglob=None, taxa=None,
                 masterTree=None,
                 output_dir=None,
                 numperm=256,
                 treeweights=None, taxfilter=None, taxmask=None,
                 lossonly=False, duplonly=False, verbose=False,
                 use_taxcodes=False,
                 datetime=datetime.now(),
                 reformat_names=False,
                 slicesubhogs=False,
                 limit_species=10, limit_events=0):
                
        """
            Initializes the LSHBuilder class with the specified parameters and sets up the necessary objects.
            
            Args:
            - tarfile_ortho (str):  path to an ensembl tarfile containing orthoxml files
            - h5_oma (str): path to an OMA hdf5 file
            - masterTree (str): path to a newick tree file
            - output_dir (str): path to the directory where the output files will be saved
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

        self.output_dir = output_dir
        self.output_dir.mkdir(parents=True, exist_ok=True)

        species = self._get_species_names()
        if masterTree is None:
            if not h5_oma:
                raise TypeError('Please specify either a database or a tree')

            self.tree_string, self.tree = phylo.get_tree(genomes=species, outdir=self.output_dir)
        else:
            # validate tree
            self.tree = phylo.TreeValidator(masterTree, species).run()

            # save the corrected tree
            phylo.to_file(self.tree, self.output_dir / 'master_tree.corrected.nwk')
            self.tree_string = phylo.to_string(self.tree)

        self.taxaIndex, self.reverse = phylo.generate_taxa_index(self.tree, self.tax_filter, self.tax_mask)
        
        with open(self.output_dir / 'taxaIndex.pkl', 'wb') as taxout:
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

        with open(self.output_dir / 'wmg.pkl', 'wb') as wmgout:
            wmgout.write( pickle.dumps(wmg))

        self.wmg = wmg

        logger.debug('\nConfiguring pyham functions')
        logger.debug('taxfilter', self.tax_filter)
        logger.debug('taxmask', self.tax_mask)
        logger.debug('swap ids', self.swap2taxcode)
        logger.debug('reformat names', self.reformat_names)
        logger.debug('use taxcodes', self.swap2taxcode)

        hamfunction = pyhamutils.get_ham_treemap_from_row
        hashfunction = hashutils.row2hash
        if slicesubhogs:
            hashfunction = hashutils.hash_trees_subhogs
        
        ### set up the pyHAM pipeline with different parameters depending on whether the input is an OMA hdf5 file or a list of orthoxml files
        logger.debug('lossonly', lossonly)
        logger.debug('duplonly', duplonly)

        self.dataset_nodes = None
        self.idmapper = None

        # remap taxfilter and taxmask 
        if taxfilter:
            #self.tax_filter = [ self.idmapper[tax] for tax in taxfilter ]
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
            #self.tax_mask = self.idmapper[taxmask]
            ### get acceptable ids here:
            tax_mask_node = self.tree.search_nodes(name=self.tax_mask)
            if tax_mask_node:
                tax_mask_node = tax_mask_node[0]
                print(f"Found tax_mask_node: {tax_mask_node.name}")
                self.dataset_nodes = [node.name for node in tax_mask_node.traverse()]
            else:
                raise RuntimeError(f"No node found with name: {self.tax_mask}")
        
        print("Kept dataset nodes:", self.dataset_nodes)

        self.HAM_PIPELINE = functools.partial(
            hamfunction, tree_string=self.tree_string, swap_ids=self.swap2taxcode,
            orthoXML_as_string=self.h5OMA is not None, slicesubhogs=self.slicesubhogs,
            limit_species=self.limit_species, limit_events=self.limit_events,
            dataset_nodes=self.dataset_nodes, verbose=self.verbose)
        ### set up the hash pipeline
        self.HASH_PIPELINE = functools.partial( hashfunction , taxaIndex=self.taxaIndex, treeweights=self.treeweights, wmg=wmg , lossonly = lossonly, duplonly = duplonly)
        print("\nSetting up input data reader")
        ### if the input is an OMA hdf5 file, set up the function to read the orthoxml data
        if self.h5OMA:
            self.READ_ORTHO = functools.partial(pyhamutils.get_orthoxml_oma, db_obj=self.db_obj)
            self.n_groups  = len(self.h5OMA.root.OrthoXML.Index)
            logger.info( 'reading oma hdf5 with n groups: %d', self.n_groups)
        elif self.fileglob:
            logger.info('reading orthoxml files: %d' , len(self.fileglob))
            self.n_groups = len(self.fileglob)
        else:
            raise RuntimeError('please specify an input file' )
        
        self.hashes_path = self.output_dir / 'hashes.h5'
        self.lshpath = self.output_dir / 'newlsh.pkl'
        self.lshforestpath = self.output_dir / 'newlshforest.pkl'
        self.mat_path = self.output_dir / 'hogmat.h5'
        self.columns = len(self.taxaIndex)
        self.verbose = verbose
        logger.debug('done')

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
                for child in oxml.iter_species(filename):
                    field = "NCBITaxId" if self.swap2taxcode else "name"
                    names.add(child.attrib[field])
                        
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
        if self.h5OMA:
            self.groups  = self.h5OMA.root.OrthoXML.Index
            ### only Fam makes sense here, the rest are OMAmer related fields
            #print(self.h5OMA.root.OrthoXML.Index.colnames)
            self.rows = len(self.groups)
            for i, row in enumerate(track_progress(
                self.groups,
                description="Processing OMA groups",
                total=self.rows,
            )):
                if i > start:
                    #### family here is HOG ID minus the "HOG:E" prefix
                    fam = row[0]
                    ### testing only
                    #if fam != 794303 and fam!= 751313: #and fam != 712236 and fam != 708323:
                    #    continue
                    ortho_fam = self.READ_ORTHO(fam)
                    ### for older versions of OMA that do not have subhogIDs
                    oxml_builder = oxml.OrthoXMLBuilder()
                    ortho_fam = oxml_builder.add_loft_ids(ortho_fam)

                    hog_size = oxml.count_species(ortho_fam)
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
            for i,file in enumerate(track_progress(
                self.fileglob,
                description="Processing OrthoXML files",
                total=len(self.fileglob),
            )):
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


    def worker(self, i, work_queue, result_queue):
        logger.debug('Starting worker #%d', i)
        while True:
            df = work_queue.get()
            if df is None:
                logger.debug("Worker #%d done", i)
                break
                
            if self.slicesubhogs:
                tree_dicts = []
                for index, row in df[["Fam", "ortho"]].iterrows():
                    result = self.HAM_PIPELINE(row)
                    tree_dicts.append(result)

                # Assign the results back to the DataFrame
                df["tree_dicts"] = tree_dicts
                df["hash_dicts"] = df[["Fam", "tree_dicts"]].apply(
                    self.HASH_PIPELINE, axis=1
                )

                # Filter out rows with empty hash_dicts to save time
                df = df[df["hash_dicts"].apply(bool)]
                if df.empty:
                    continue

                result_dict = {}
                for _, row in df.iterrows():
                    for subhog in row["hash_dicts"]:
                        # don't save orthoxml strings if OMA
                        ortho = row["ortho"] if self.fileglob else ""

                        result_dict[(row["Fam"], subhog)] = {
                            "tree": row["tree_dicts"][subhog],
                            "hash": row["hash_dicts"][subhog][1],
                            "ortho": ortho
                        }
                result_dict = pd.DataFrame.from_dict(result_dict, orient="index")
                result_queue.put(result_dict)
            else:
                df["tree"] = df[["Fam", "ortho"]].apply(self.HAM_PIPELINE, axis=1)
                # add a dictionary of results with subhogs { fam_sub1: { 'tree':tp , 'Fam':fam }  , fam_sub2: { 'tree':tp , 'Fam':fam } , ... }
                # returned_df = pd.DataFrame.from_dict(df['tree'].to_dict(), orient='index')
                # merge with pandas on right e.g. df.merge( returned_df , on = 'Fam' , how = 'right' )

                df[["hash", "rows"]] = df[["Fam", "tree"]].apply(self.HASH_PIPELINE, axis=1)
                if self.fileglob:
                    result_queue.put(df[["Fam", "hash", "ortho"]])
                else:
                    result_queue.put(df[["Fam", "hash"]])

    
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
        mapping.to_csv(self.output_dir / 'fam2orthoxml.csv')


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

        total_subfam_ids = []
        totals_subfams = 0

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

                        # handle slicesubhogs
                        if self.slicesubhogs:
                            subfam_ids_list = this_dataframe.index.to_list()
                            total_subfam_ids.extend(subfam_ids_list)

                            nsubfams = len(hashes)
                            # Resize if necessary
                            if h5hashes[taxstr].shape[0] < totals_subfams + nsubfams + chunk_size:
                                example = next(iter(hashes.keys()))
                                h5hashes[taxstr].resize((totals_subfams + nsubfams + chunk_size,
                                                         len(hashes[example].hashvalues.ravel())))

                            # Store hashes
                            h5hashes[taxstr][totals_subfams:totals_subfams + nsubfams, :] = [hashes[fam].hashvalues.ravel() for fam in hashes]

                            # Update forest
                            for fam in hashes:
                                forest.add(str(fam[0]) + '_' + str(fam[1]), hashes[fam])

                            totals_subfams += nsubfams
                        else:
                            hashes = { fam:hashes[fam] for fam in hashes if hashes[fam]}
                            for fam in hashes:
                                forest.add(str(fam), hashes[fam])

                            for fam in hashes:
                                # if all rows are filled, allocate the next chunk
                                if dataset.shape[0] <= fam:
                                    dataset.resize(fam + chunk_size, axis=0)

                                dataset[fam, :] = hashes[fam].hashvalues.ravel()

                        count += len(hashes)
                        if self.fileglob:
                            mapping_frames.append(this_dataframe[['ortho']])

                        if hashes and t.time() - save_start > 200:
                            h5hashes.flush()
                            self._index_and_save(forest)
                            logger.debug('Testing forest')
                            logger.debug(forest.query(hashes[fam], k=10))

                            if self.fileglob:
                                self._save_family_mapping(mapping_frames)
                            save_start = t.time()
                    else:
                        logger.debug('Saver received an empty result batch')
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
                csv_path = os.path.join(self.output_dir, 'fam2orthoxml.csv')
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
    
    def run_pipeline(self , threads):
        logger.info('Running with %d threads', threads)
        # Levels mode counts represented species on each parsed HOG
        minhog_size = None if self.slicesubhogs else self.limit_species
        data_generator = self.generates_dataframes(size=100, minhog_size=minhog_size)
        self._run_parallel(threads=threads, data_generator=data_generator)
        return self.hashes_path, self.lshforestpath, self.mat_path


def run(parsed_args):
    """Build profiles from arguments validated by the lightweight CLI."""

    taxfilter = None
    taxmask = None
    omafile = None

    orthoglob = None
    args = vars(parsed_args)

    if 'OrthoGlob' in args:
        if args['OrthoGlob']:
            print("Using orthoxml files from:", args['OrthoGlob'])
            orthoglob = glob.glob(args['OrthoGlob'])
    
    if args['dbtype']:
        taxfilter = DB_PRESETS[args['dbtype']]['taxfilter']
        taxmask = DB_PRESETS[args['dbtype']]['taxmask']
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

    _args = parsed_args
    output_dir = _args.outpath
    njobs = _args.njobs
    if _args.nthreads:
        logger.warning("--nthreads is deprecated. Please use --njobs")
        njobs = _args.nthreads

    if args['taxweights']:
        from keras.models import model_from_json
        json_file = open(  args['taxweights']+ '.json', 'r')
        loaded_model_json = json_file.read()
        json_file.close()
        model = model_from_json(loaded_model_json)
        # load weights into new model
        model.load_weights(  args['taxweights']+".h5")
        logger.info("Loaded model from disk")
        weights = model.get_weights()[0]
        weights += 10 ** -10
    else:
        weights = None
    if args['mastertree']:
        mastertree = Path(args['mastertree'])
    else:
        mastertree=None

    start = time.time()
    if omafile:
        with open_file(omafile , mode="r") as h5_oma:
            lsh_builder = LSHBuilder(h5_oma=h5_oma, fileglob=orthoglob, output_dir=output_dir,
                                     numperm=nperm,
                                     treeweights= weights , taxfilter = taxfilter, taxmask=taxmask , masterTree =mastertree ,
                                     lossonly=_args.lossonly , duplonly = _args.duplonly , use_taxcodes = _args.taxcodes , reformat_names=_args.reformat_names,
                                     slicesubhogs=_args.slicesubhogs, limit_species=args['specieslim'], limit_events=args['eventslim'],
                                     verbose=_args.verbose,)
            lsh_builder.run_pipeline(njobs)

    else:
        lsh_builder = LSHBuilder(h5_oma = None,  fileglob=orthoglob, output_dir=output_dir,
                                 numperm = nperm,
                                 treeweights= weights , taxfilter = taxfilter, taxmask=taxmask ,
                                 masterTree =mastertree , lossonly = _args.lossonly , duplonly = _args.duplonly , use_taxcodes = _args.taxcodes ,
                                 reformat_names=_args.reformat_names,
                                 slicesubhogs=_args.slicesubhogs, limit_species=args['specieslim'],
                                 limit_events=args['eventslim'],
                                 verbose=_args.verbose)
        lsh_builder.run_pipeline(njobs)

    logger.info("Done in %.2fs", time.time() - start)