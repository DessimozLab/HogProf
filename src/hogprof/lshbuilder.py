import functools
import glob
import logging
import multiprocessing as mp
import pickle
import queue
import random
import time
import time as t
import xml.etree.ElementTree as ET
from contextlib import nullcontext
from datetime import datetime
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import polars as pl
from datasketch import MinHashLSHForest, WeightedMinHash, WeightedMinHashGenerator
from pyoma.browser import db
from tables import open_file

from hogprof.artifacts import FOREST, HASHSTORE, MAPPING, MINHASH_GENERATOR, TAXA_INDEX, TREE
from hogprof import hashstore
from hogprof.cli import track_progress
from hogprof.hashstore import HashStore
from hogprof.manifest import IndexManifest
from hogprof.metadata import HogEntry, HogEntryIndex
from hogprof.presets import DB_PRESETS
from hogprof.utils import hashutils, phylo, pyhamutils
from hogprof.utils import orthoxml as oxml

logger = logging.getLogger(__name__)


class LSHBuilder:
    """
    Makes a phylogenetic profiling database from input orthoxml files
    a taxonomic tree. Takes an OMA hdf5 file or a glob expression
    containing orthoxml file with orthologous groups as input.
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
                 limit_species=10, limit_events=0,
                 seed=0, provenance=None):
                
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
            - limit_events (int): the minimum number of events (losses + duplications) a profile must have to be included in the database; 0 keeps all (default: 0)
            - provenance (dict): extra information about the inputs to record in the index manifest

        """
        logger.debug("Initializing LSHBuilder")

        random.seed(seed)
        np.random.seed(seed)

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
        self.lossonly = lossonly
        self.duplonly = duplonly
        self.seed = seed
        self.provenance = dict(provenance or {})

        self.min_hog_size = 0 if self.slicesubhogs else self.limit_species
        self.max_hog_size = np.inf

        self.output_dir = output_dir
        self.output_dir.mkdir(parents=True, exist_ok=True)
        # an old manifest would mark a partially rebuilt index as complete
        IndexManifest.remove(self.output_dir)

        species = self._get_species_names()
        self.species_rename = {}
        if masterTree is not None and self.h5OMA and not self.swap2taxcode:
            self.species_rename = self._oma_species_rename(masterTree)
            species = {self.species_rename.get(name, name) for name in species}

        if masterTree is None:
            if not h5_oma:
                raise TypeError('Please specify either a database or a tree')

            self.tree_string, self.tree = phylo.get_tree(genomes=species, outdir=self.output_dir)
        else:
            # validate tree
            self.tree = phylo.TreeValidator(masterTree, species).run()

            # save the corrected tree
            phylo.to_file(self.tree, self.output_dir / TREE)
            self.tree_string = phylo.to_string(self.tree)

        self.taxaIndex, self.reverse = phylo.generate_taxa_index(self.tree, self.tax_filter, self.tax_mask)
        
        with open(self.output_dir / TAXA_INDEX, 'wb') as taxout:
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
        wmg = WeightedMinHashGenerator(3*tax_max , sample_size = numperm , seed=hashutils.HASH_SEED)

        with open(self.output_dir / MINHASH_GENERATOR, 'wb') as wmgout:
            wmgout.write( pickle.dumps(wmg))

        self.wmg = wmg

        logger.debug("Configuring pyham functions")
        logger.debug("taxfilter: %s", self.tax_filter)
        logger.debug("taxmask: %s", self.tax_mask)
        logger.debug("swap ids: %s", self.swap2taxcode)
        logger.debug("reformat names: %s", self.reformat_names)
        logger.debug("use taxcodes: %s", self.swap2taxcode)

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

        self.HAM_PIPELINE = functools.partial(
            hamfunction, tree_string=self.tree_string, swap_ids=self.swap2taxcode,
            orthoXML_as_string=self.h5OMA is not None, slicesubhogs=self.slicesubhogs,
            limit_species=self.limit_species, limit_events=self.limit_events,
            dataset_nodes=self.dataset_nodes, verbose=self.verbose)

        ### set up the hash pipeline
        self.HASH_PIPELINE = functools.partial( hashfunction , taxaIndex=self.taxaIndex, treeweights=self.treeweights, wmg=wmg , lossonly = lossonly, duplonly = duplonly)

        logger.debug("Setting up input data reader")
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
        
        self.hashes_path = self.output_dir / HASHSTORE
        self.lshpath = self.output_dir / 'newlsh.pkl'
        self.lshforestpath = self.output_dir / FOREST
        self.mat_path = self.output_dir / 'hogmat.h5'
        self.columns = len(self.taxaIndex)
        self.verbose = verbose
        logger.debug('Done.')

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

    def _oma_species_rename(self, tree):
        """
        Some OMA species trees (e.g. All.Jul2024) uses scientific
        names of the species in both a leaf and an internal node.
        Map those internal nodes to their Uniprot codes
        """
        tree = phylo.from_file(tree)
        leaves = set(tree.get_leaf_names())
        internal = {node.name for node in tree.traverse() if not node.is_leaf()}
        genomes = self.h5OMA.root.Genome.read()

        rename = {}
        for name, code in zip(genomes["SciName"], genomes["UniProtSpeciesCode"]):
            name, code = name.decode(), code.decode()
            if name in internal and code in leaves:
                rename[name] = code

        if rename:
            logger.warning(
                "Using UniProt codes for %d species whose names "
                "are internal tree nodes: %s",
                len(rename),
                ", ".join(f"{k} -> {v}" for k, v in sorted(rename.items())),
            )
        return rename

    def _rename_species(self, xml):
        """Apply species_rename to the <species> elements of a parsed OrthoXML"""
        if self.species_rename:
            for node in xml.findall("{http://orthoXML.org/2011/}species"):
                name = node.get("name")
                if name in self.species_rename:
                    node.set("name", self.species_rename[name])
        return xml

    def load_one(self, fam):
        #test function to try out the pipeline on one orthoxml
        ortho_fam = self.READ_ORTHO(fam)
        pyham_tree = self.HAM_PIPELINE([fam, ortho_fam])
        hog_matrix,weighted_hash = hashutils.hash_tree(pyham_tree , self.taxaIndex , self.treeweights , self.wmg)
        return ortho_fam , pyham_tree, weighted_hash,hog_matrix

    def _hog_size_ok(self, hog_size):
        result = (self.max_hog_size is None or hog_size <= self.max_hog_size)
        result &= (self.min_hog_size is None or hog_size >= self.min_hog_size)
        return result

    def _check_hog_size(self, orthoxml_file):
        with open(orthoxml_file) as ortho:
            ortho_str = ortho.read()

        # filtering by HOG size in the # of species
        hog_size = oxml.count_species(ortho_str)
        return self._hog_size_ok(hog_size)

    def generate_dataframes(self, size=100):
        families = {}
        start = -1
        if self.h5OMA:
            self.groups  = self.h5OMA.root.OrthoXML.Index
            buffer = self.h5OMA.root.OrthoXML.Buffer
            oxml_builder = oxml.OrthoXMLBuilder()
            self.rows = len(self.groups)

            unique_hog_id = 0

            for i, row in enumerate(track_progress(
                self.groups,
                description="Processing OMA groups",
                total=self.rows,
            )):
                if i > start:
                    # The current index row already locates the unaugmented XML
                    offset = int(row['HogBufferOffset'])
                    stop = offset + int(row['HogBufferLength'])
                    xml = ET.fromstring(buffer[offset:stop].tobytes().decode())

                    # filtering by HOG size in the # of species
                    hog_size = oxml.count_species(xml)
                    if self._hog_size_ok(hog_size):
                        ortho_fam = oxml_builder.add_loft_ids(self._rename_species(xml))
                        families[unique_hog_id] = {
                            'ortho': ortho_fam,
                            'hog_id': oxml.get_root_hog_id(xml),
                        }
                        unique_hog_id += 1

                    if len(families) >= size:
                        pd_dataframe = pd.DataFrame.from_dict(families, orient='index')
                        pd_dataframe['hog_key'] = pd_dataframe.index
                        yield pd_dataframe
                        families = {}
            if families:
                pd_dataframe = pd.DataFrame.from_dict(families, orient='index')
                pd_dataframe['hog_key'] = pd_dataframe.index
                yield pd_dataframe

        elif self.fileglob:
            unique_fam_id = 0

            for orthoxml_file in track_progress(
                self.fileglob,
                description="Processing OrthoXML files",
                total=len(self.fileglob),
            ):
                if self._check_hog_size(orthoxml_file):
                    families[unique_fam_id] = {'ortho': orthoxml_file, 'hog_id': ''}
                    unique_fam_id += 1

                if len(families) >= size:
                    pd_dataframe = pd.DataFrame.from_dict(families, orient='index')
                    pd_dataframe['hog_key'] = pd_dataframe.index
                    #print(pd_dataframe)
                    yield pd_dataframe
                    families = {}

            if families:
                pd_dataframe = pd.DataFrame.from_dict(families, orient='index')
                pd_dataframe['hog_key'] = pd_dataframe.index
                yield pd_dataframe

    def _extract_hog_name(self, orthoxml_file: Path) -> str:
        return Path(orthoxml_file).name.split(".", 1)[0]

    def worker(self, i, work_queue, result_queue):
        logger.debug('Starting worker #%d', i)
        while True:
            df = work_queue.get()
            if df is None:
                logger.debug("Worker #%d done", i)
                break
                
            if self.slicesubhogs:
                # Hash every subHOG profile as soon as it's made and drop it.
                # Minimizes peak RAM consumption
                rows = []
                for _, row in df.iterrows():
                    fam_id = self._extract_hog_name(row["ortho"]) if self.fileglob else str(row["hog_id"])
                    tree_maps = self.HAM_PIPELINE(
                        (row["hog_key"], row["ortho"]), lazy=True
                    )
                    for subhog, tree_map in tree_maps:
                        level = tree_map.name
                        _, hashed = hashutils.hash_tree(
                            tree_map,
                            self.taxaIndex,
                            self.treeweights,
                            self.wmg,
                            lossonly=self.lossonly,
                            duplonly=self.duplonly,
                        )
                        if hashed is None:
                            continue

                        rows.append(
                            {
                                "hog_key": row["hog_key"],
                                "fam_id": fam_id,
                                "hog_id": subhog.removeprefix(f"{level}_"),
                                "level": level,
                                # Polars Object columns cannot be pickled by mp.Queue.
                                # Send native values and reconstruct the hash in saver
                                "hash": hashed.hashvalues.tolist(),
                                "hash_seed": hashed.seed,
                                "source": str(row["ortho"]) if self.fileglob else "",
                            }
                        )
                if rows:
                    result_queue.put(pl.DataFrame(rows))
            else:
                df["tree"] = df[["hog_key", "ortho"]].apply(self.HAM_PIPELINE, axis=1)
                df[["hash", "rows"]] = df[["hog_key", "tree"]].apply(self.HASH_PIPELINE, axis=1)

                # HOG names are extracted from the names of orthoxml files.
                # Let's hope they are correct
                if self.fileglob:
                    df["hog_id"] = df["ortho"].apply(self._extract_hog_name)
                    df["source"] = df["ortho"].map(str)
                else:
                    df["hog_id"] = df["hog_id"].map(str)
                    df["source"] = ""
                df["fam_id"] = df["hog_id"]
                df["level"] = df["tree"].map(lambda tree: tree.name if tree is not None else "")

                result_queue.put(df[[*HogEntryIndex.COLUMNS, "hash"]])

    
    def _make_tax_str(self):
        taxstr = ""
        if self.tax_filter is None:
            taxstr = "NoFilter"
        if self.tax_mask is None:
            taxstr += "NoMask"
        else:
            taxstr = str(self.tax_filter)
        return taxstr

    def _save_forest(self, forest):
        logger.debug("Saving forest")

        forest.index()
        with open(self.lshforestpath, "wb") as forest_out:
            forest_out.write(pickle.dumps(forest, -1))

        logger.debug("Forest saved")

    def _forest_from_hashes(self, block_size=65536):
        """Build the LSH forest from the saved hashes, keyed by their rows"""
        forest = MinHashLSHForest(num_perm=self.numperm)
        with HashStore(self.hashes_path, num_permutations=self.numperm,
                       dataset=self._make_tax_str(), seed=self.wmg.seed) as store:
            for start in range(0, store.num_rows, block_size):
                block = store.get(np.arange(start, min(start + block_size, store.num_rows)))
                for offset, packed in enumerate(block):
                    forest.add(str(start + offset), hashstore.to_minhash(packed, seed=store.seed))
        return forest

    def _save_metadata(self, entries):
        df = pd.DataFrame(entries, columns=HogEntryIndex.COLUMNS)
        df.to_csv(self.output_dir / MAPPING, index=False)


    def saver(self, result_queue):
        global_time = t.time()
        self.start_time = global_time

        chunk_size = 100
        count = 0
        # With --slicesubhogs, the index is reordered by level after all
        # hashes are saved, which changes hog_keys: the forest is built
        # from the final hashes then (see `run_pipeline`).
        forest = None if self.slicesubhogs else MinHashLSHForest(num_perm=self.numperm)
        metadata_entries = []
        taxstr = self._make_tax_str()

        # With --slicesubhogs, we don't know in advance into how many
        # subHOGs a family will be expanded. Therefore, the global hog_key
        # has to be reassigned so that:
        #   1) every subHOG gets a unique key
        #   2) hog_keys are continuous
        hog_key = 0

        with h5py.File(self.hashes_path, 'w', libver='latest') as h5hashes:
            hash_width = 2 * self.numperm

            # The main hash mapping. LSHForest.add takes a pair = (key, MihHash object),
            # however you cannot query it with 'key', only with MinHash.
            # To support that, we store the mapping { key -> MinHash }
            # for everything that is indexed in the forest.
            dataset = h5hashes.create_dataset(
                taxstr,
                (chunk_size, hash_width),
                maxshape=(None, hash_width),
                dtype="int32",
            )
            h5hashes.flush()
            logger.debug('Creating dataset filtered at taxonomic level: %s', taxstr)

            while True:
                current_df = result_queue.get()
                if current_df is None:
                    break

                # Temporary: Pandas is to be replaced by Polars everywhere.
                # currently we still use pd for the rootHOG mode
                is_polars = isinstance(current_df, pl.DataFrame)
                rows = current_df.to_dicts() if is_polars else current_df.to_dict("records")
                for row in rows:
                    hashed = row["hash"]
                    if hashed is None:
                        continue

                    if is_polars:
                        hashed = WeightedMinHash(
                            seed=row["hash_seed"],
                            hashvalues=np.asarray(hashed, dtype=np.int64),
                        )

                    key = hog_key if self.slicesubhogs else int(row["hog_key"])

                    # resize the HogEntryIndex dataset if all rows filled
                    if dataset.shape[0] <= key:
                        dataset.resize(key + chunk_size, axis=0)

                    # save the hash
                    dataset[key, :] = hashed.hashvalues.ravel()

                    # index the hash
                    if forest is not None:
                        forest.add(str(key), hashed)

                    metadata_row = HogEntry(
                        hog_key=key, fam_id=row["fam_id"], hog_id=row["hog_id"],
                        level=row["level"], source=row["source"],
                    )
                    metadata_entries.append(metadata_row)
                    hog_key += 1
                    count += 1

            logger.info('Wrapping up the run')
            if forest is not None:
                self._save_forest(forest)
            h5hashes.flush()
            self._save_metadata(metadata_entries)
            logger.info('Saver wrote %d hashes', count)


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


    def _save_manifest(self):
        num_hogs = (
            pl.scan_csv(self.output_dir / MAPPING)
            .select(pl.len()).collect().item()
        )
        manifest = IndexManifest(
            num_permutations=self.numperm,
            hash_seed=self.wmg.seed,
            hashes_dataset=self._make_tax_str(),
            slice_subhogs=self.slicesubhogs,
            num_hogs=num_hogs,
            build={
                **self.provenance,
                "taxon_mask": self.tax_mask,
                "exclude_taxa": self.tax_filter,
                "min_species": self.limit_species,
                "min_events": self.limit_events,
                "loss_only": self.lossonly,
                "duplication_only": self.duplonly,
                "use_tax_ids": self.swap2taxcode,
                "seed": self.seed,
            },
        )
        path = manifest.write(self.output_dir)
        logger.info('Wrote index manifest: %s', path)

    def run_pipeline(self , threads):
        logger.info('Running with %d threads', threads)
        # Levels mode counts represented species on each parsed HOG
        data_generator = self.generate_dataframes(size=100)
        self._run_parallel(threads=threads, data_generator=data_generator)

        if self.slicesubhogs:
            # Hashes were computed family by family, so every level is
            # scattered over the whole file. Make levels contiguous, so that
            # level-wide searches read one block.
            logger.info('Reordering the index by level')
            hashstore.reorder_by_level(self.hashes_path, self.output_dir / MAPPING,
                                       dataset=self._make_tax_str())
            logger.info('Building the LSH forest')
            self._save_forest(self._forest_from_hashes())

        # written last: its presence marks a complete index
        self._save_manifest()
        return self.hashes_path, self.lshforestpath, self.mat_path


def run(
    *,
    output_dir,
    tree=None,
    oma=None,
    orthoxml_glob=None,
    tax_weights=None,
    taxon_mask=None,
    exclude_taxa=None,
    min_species=10,
    min_events=0,
    db_type=None,
    num_permutations=256,
    loss_only=False,
    duplication_only=False,
    use_tax_ids=False,
    slice_subhogs=False,
    seed=0,
    njobs=1,
    verbose=False,
):
    if (oma is not None) == (orthoxml_glob is not None):
        raise ValueError("Supply exactly one of oma and orthoxml_glob")

    if loss_only and duplication_only:
        raise ValueError("loss_only and duplication_only are mutually exclusive")

    if db_type is not None:
        preset = DB_PRESETS[db_type]
        if taxon_mask is None:
            taxon_mask = preset['taxon_mask']
        if exclude_taxa is None:
            exclude_taxa = preset['exclude_taxa']

    weights = None
    if tax_weights is not None:
        from keras.models import model_from_json

        with open(f"{tax_weights}.json") as model_file:
            model = model_from_json(model_file.read())
        model.load_weights(f"{tax_weights}.h5")
        logger.info("Loaded model from disk")
        weights = model.get_weights()[0] + 1e-10

    # Discover input files from the glob
    orthoxml_files = None
    if orthoxml_glob is not None:
        orthoxml_files = glob.glob(orthoxml_glob)
        logger.info("Using OrthoXML files from: %s", orthoxml_glob)

    start = time.monotonic()
    database = open_file(oma, mode="r") if oma is not None else nullcontext(None)
    with database as h5_oma:
        lsh_builder = LSHBuilder(
            h5_oma=h5_oma,
            fileglob=orthoxml_files,
            output_dir=Path(output_dir),
            masterTree=Path(tree) if tree is not None else None,
            numperm=num_permutations,
            treeweights=weights,
            taxfilter=exclude_taxa,
            taxmask=taxon_mask,
            lossonly=loss_only,
            duplonly=duplication_only,
            use_taxcodes=use_tax_ids,
            slicesubhogs=slice_subhogs,
            limit_species=min_species,
            limit_events=min_events,
            seed=seed,
            verbose=verbose,
            provenance={
                "oma": str(oma) if oma is not None else None,
                "orthoxml_glob": orthoxml_glob,
                "species_tree": str(tree) if tree is not None else None,
                "db_type": db_type,
                "tax_weights": tax_weights,
            },
        )
        paths = lsh_builder.run_pipeline(njobs)

    logger.info("Done in %.2fs", time.monotonic() - start)
    return paths
