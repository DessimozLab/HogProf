import argparse
import functools
import gc
import logging
import multiprocessing as mp
import os
import pickle
import random
import time
import time as t
from datetime import datetime

import ete3
import h5py
import numpy as np
import pandas as pd
import tqdm
from datasketch import MinHashLSHForest, WeightedMinHashGenerator
from pyoma.browser import db
from tables import *

from HogProf.utils import pyhamutils, hashutils, files_utils

random.seed(0)
np.random.seed(0)

logging.basicConfig(level=logging.DEBUG,
                    format='%(asctime)x [%(levelname)s] %(message)s',
                    handlers=[
                        logging.FileHandler('debug.log'),
                        logging.StreamHandler()])


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

    def __init__(self, h5_oma=None, fileglob=None, taxa=None, masterTree=None, saving_name=None, numperm=256,
                 treeweights=None, taxfilter=None, taxmask=None, lossonly=False, duplonly=False, verbose=False,
                 use_taxcodes=False, datetime=datetime.now()):

        """
            Initializes the LSHBuilder class with the specified parameters and sets up the necessary objects.

            Args:
            - tarfile_ortho (str):  path to an ensembl tarfile containing orthoxml files
            - h5_oma (str): path to an OMA hdf5 file
            - taxa (str): path to a file containing a list of taxonomic codes for all the species in the db
            - masterTree (str): path to a newick tree file
            - saving_name (str): path to the directory where the output files will be saved
            - numperm (int): the number of permutations to use in the MinHash generation (default: 256)
            - treeweights (str): path to a pickled file containing the weights for the tree
            - taxfilter (str): path to a file containing a list of taxonomic codes to filter from the tree
            - taxmask (str): path to a file containing a list of taxonomic codes to mask from the tree
            - verbose (bool): whether to print verbose output (default: False)

        """
        logging.info('Initialising %s' % self.__class__.__name__)
        self.groups = None
        self.errorfile = None

        if h5_oma:
            self.h5OMA = h5_oma
            self.db_obj = db.Database(h5_oma)
            self.oma_id_obj = db.OmaIdMapper(self.db_obj)
        else:
            self.h5OMA = None
            self.db_obj = None
            self.oma_id_obj = None

        self.tax_filter = taxfilter
        self.tax_mask = taxmask
        self.verbose = verbose
        self.datetime = datetime
        self.fileglob = fileglob
        self.date_string = "{:%B_%d_%Y_%H_%M}".format(datetime.now())

        if saving_name:
            self.saving_name = saving_name
            if self.saving_name[-1] != '/':
                self.saving_name = self.saving_name + '/'
            self.saving_path = saving_name
            if not os.path.isdir(self.saving_path):
                os.mkdir(path=self.saving_path)
        else:
            raise Exception('please specify an output location')

        if masterTree is None:
            if h5_oma:
                genomes = pd.DataFrame(h5_oma.root.Genome.read())["NCBITaxonId"].tolist()
                genomes = [str(g) for g in genomes]
                taxa = genomes + [131567, 2759, 2157, 45596] + [taxrel[0] for taxrel in
                                                                list(h5_oma.root.Taxonomy[:])] + [taxrel[1] for taxrel
                                                                                                  in list(
                        h5_oma.root.Taxonomy[:])]
                self.tree_string, self.tree_ete3 = files_utils.get_tree(taxa=taxa, genomes=genomes,
                                                                        outdir=self.saving_path)
            else:
                raise Exception('please specify either a list of taxa or a tree')
            self.swap2taxcode = True
        elif masterTree:
            self.tree_ete3 = ete3.Tree(masterTree, format=1)
            with open(masterTree) as treein:
                self.tree_string = treein.read()
            self.swap2taxcode = use_taxcodes

        self.taxaIndex, self.reverse = files_utils.generate_taxa_index(self.tree_ete3, self.tax_filter, self.tax_mask)

        with open(self.saving_path + 'taxaIndex.pkl', 'wb') as taxout:
            taxout.write(pickle.dumps(self.taxaIndex))

        self.numperm = numperm

        if treeweights is None:
            # generate aconfig_utilsll ones
            self.treeweights = hashutils.generate_treeweights(self.tree_ete3, self.taxaIndex, taxfilter, taxmask)
        else:
            # load machine learning weights
            self.treeweights = treeweights

        wmg = WeightedMinHashGenerator(3 * len(self.taxaIndex), sample_size=numperm, seed=1)

        with open(self.saving_path + 'wmg.pkl', 'wb') as wmgout:
            wmgout.write(pickle.dumps(wmg))
        self.wmg = wmg

        logging.info('Configuring pyham functions')

        if self.h5OMA:
            self.HAM_PIPELINE = functools.partial(pyhamutils.get_ham_treemap_from_row, tree=self.tree_string,
                                                  swap_ids=self.swap2taxcode)
        else:
            self.HAM_PIPELINE = functools.partial(pyhamutils.get_ham_treemap_from_row, tree=self.tree_string,
                                                  swap_ids=self.swap2taxcode, orthoXML_as_string=False)

        self.HASH_PIPELINE = functools.partial(hashutils.row2hash, taxaIndex=self.taxaIndex,
                                               treeweights=self.treeweights, wmg=wmg, lossonly=lossonly,
                                               duplonly=duplonly)

        if self.h5OMA:
            self.READ_ORTHO = functools.partial(pyhamutils.get_orthoxml_oma, db_obj=self.db_obj)

        if self.h5OMA:
            self.n_groups = len(self.h5OMA.root.OrthoXML.Index)
        elif self.fileglob:
            self.n_groups = len(self.fileglob)
        else:
            raise Exception('please specify an input file')

        self.hashes_path = self.saving_path + 'hashes.h5'
        self.lshpath = self.saving_path + 'newlsh.pkl'
        self.lshforestpath = self.saving_path + 'newlshforest.pkl'
        self.mat_path = self.saving_path + 'hogmat.h5'
        self.columns = len(self.taxaIndex)
        self.verbose = verbose
        print('Initialised')

    def load_one(self, fam):
        # test function to try out the pipeline on one orthoxml
        ortho_fam = self.READ_ORTHO(fam)
        pyham_tree = self.HAM_PIPELINE([fam, ortho_fam])
        hog_matrix, weighted_hash = hashutils.hash_tree(pyham_tree, self.taxaIndex, self.treeweights, self.wmg)
        return ortho_fam, pyham_tree, weighted_hash, hog_matrix

    def generates_dataframes(self, size=100, minhog_size=10, maxhog_size=None):
        families = {}
        start = -1
        if self.h5OMA:
            self.groups = self.h5OMA.root.OrthoXML.Index
            self.rows = len(self.groups)

            for i, row in tqdm.tqdm(enumerate(self.groups)):
                if i > start:
                    fam = row[0]
                    ortho_fam = self.READ_ORTHO(fam)
                    hog_size = ortho_fam.count('<species name=')

                    if (maxhog_size is None or hog_size < maxhog_size) and (
                            minhog_size is None or hog_size > minhog_size):
                        families[fam] = {'ortho': ortho_fam}

                    if len(families) > size:
                        pd_dataframe = pd.DataFrame.from_dict(families, orient='index')
                        pd_dataframe['Fam'] = pd_dataframe.index
                        yield pd_dataframe
                        families = {}

            pd_dataframe = pd.DataFrame.from_dict(families, orient='index')
            pd_dataframe['Fam'] = pd_dataframe.index

            yield pd_dataframe

            logging.info('Last dataframe sent')

        elif self.fileglob:
            for i, file in enumerate(tqdm.tqdm(self.fileglob)):

                with open(file) as ortho:
                    orthostr = ortho.read()

                hog_size = orthostr.count('<species name=')

                if (maxhog_size is None or hog_size < maxhog_size) and (minhog_size is None or hog_size > minhog_size):
                    families[i] = {'ortho': file}

                if len(families) > size:
                    pd_dataframe = pd.DataFrame.from_dict(families, orient='index')
                    pd_dataframe['Fam'] = pd_dataframe.index
                    yield pd_dataframe
                    families = {}

                if i % 10000 == 0:
                    print(i)
                    # save the mapping of fam to orthoxml
                    pd_dataframe = pd.DataFrame.from_dict(families, orient='index')
                    pd_dataframe['Fam'] = pd_dataframe.index
                    pd_dataframe.to_csv(self.saving_path + 'fam2orthoxml.csv')

    def worker(self, i, q, retq, matq, l):
        if self.verbose:
            logging.info('Initialising worker %s ' % str(i))
        while True:
            df = q.get()
            if df is not None:
                df['tree'] = df[['Fam', 'ortho']].apply(self.HAM_PIPELINE, axis=1)
                df[['hash', 'rows']] = df[['Fam', 'tree']].apply(self.HASH_PIPELINE, axis=1)
                retq.put(df[['Fam', 'hash']])
            else:
                if self.verbose:
                    print('Worker done %s' % str(i))
                break

    def saver(self, i, q, retq, matq, l):
        save_start = t.time()
        global_time = t.time()
        chunk_size = 100
        count = 0
        forest = MinHashLSHForest(num_perm=self.numperm)
        taxstr = ''

        if self.tax_filter is None:
            taxstr = 'NoFilter'

        if self.tax_mask is None:
            taxstr += 'NoMask'
        else:
            taxstr = str(self.tax_filter)

        self.errorfile = self.saving_path + 'errors.txt'
        with open(self.errorfile, 'w') as hashes_error_files:
            with h5py.File(self.hashes_path, 'w', libver='latest') as h5hashes:
                datasets = {}

                if taxstr not in h5hashes.keys():
                    if self.verbose:
                        logging.info('Creating dataset')
                        logging.info('Filtered at taxonomic level: ' + taxstr)
                    h5hashes.create_dataset(taxstr, (chunk_size, 0), maxshape=(None, None), dtype='int32')

                    if self.verbose:
                        logging.info(datasets)
                    h5flush = h5hashes.flush

                logging.info('Initialising saver %s ' % str(i))

                while True:
                    this_dataframe = retq.get()
                    if this_dataframe is not None:
                        if not this_dataframe.empty:
                            hashes = this_dataframe['hash'].to_dict()
                            hashes = {fam: hashes[fam] for fam in hashes if hashes[fam]}
                            [forest.add(str(fam), hashes[fam]) for fam in hashes]

                            for fam in hashes:
                                if len(h5hashes[taxstr]) < fam + 10:
                                    h5hashes[taxstr].resize((fam + chunk_size, len(hashes[fam].hashvalues.ravel())))
                                h5hashes[taxstr][fam, :] = hashes[fam].hashvalues.ravel()
                                count += 1

                            if t.time() - save_start > 200:
                                logging.info(t.time() - global_time)
                                forest.index()
                                logging.info(forest.query(hashes[fam], k=10))
                                h5flush()
                                save_start = t.time()

                                with open(self.lshforestpath, 'wb') as forestout:
                                    forestout.write(pickle.dumps(forest, -1))

                                if self.verbose:
                                    logging.info('Saved to %s' % str(t.time() - global_time))
                        else:
                            print(this_dataframe)
                    else:
                        if self.verbose:
                            logging.info('Wrapping it up')

                        with open(self.lshforestpath, 'wb') as forestout:
                            forestout.write(pickle.dumps(forest, -1))

                        h5flush()

                        if self.verbose:
                            print('DONE SAVER' + str(i))
                        break

    def matrix_updater(self, iprocess, q, retq, matq, l):
        logging.info('Initialising hogmat saver ' + str(iprocess))
        h5mat = None
        times1 = []
        frames = []
        with h5py.File(self.mat_path + str(iprocess) + 'h5', 'w', libver='latest') as h5hashes:
            i = 0
            while True:
                rows = matq.get()
                if rows is not None:
                    rows = rows.dropna()
                    maxfam = rows.Fam.max()
                    if h5mat is None:
                        h5hashes.create_dataset('matrows', (10, block.shape[1]), maxshape=(None, block.shape[1]),
                                                chunks=(1, block.shape[1]), dtype='i8')
                        h5mat = h5hashes['matrows']
                    if h5mat.shape[0] < maxfam:
                        h5mat.resize((maxfam + 1, block.shape[1]))
                    i += 1
                    frames.append(rows)
                    assign = t.time()
                    index = np.asarray(rows.Fam)
                    block = np.vstack(rows.rows)
                    h5mat[index, :] = block

                    times1.append(t.time() - assign)
                    if len(times1) > 10:
                        times1.pop(0)
                        logging.info('Mean time: %s' % np.mean(times1))
                    h5hashes.flush()
                else:
                    h5hashes.flush()
                    break
        logging.info('DONE MAT UPDATER %s' % str(i))

    def run_pipeline(self, threads):
        logging.info('Running with %s threads:' % threads)
        functype_dict = {'worker': (self.worker, threads, True), 'updater': (self.saver, 1, False),
                         'matrix_updater': (self.matrix_updater, 0, False)}

        def mp_with_timeout(functypes, data_generator):
            lock = mp.Lock()
            cores = mp.cpu_count()
            q = mp.Queue(maxsize=cores * 10)
            retq = mp.Queue(maxsize=cores * 10)
            matq = mp.Queue(maxsize=cores * 10)
            work_processes = {}
            logging.info('Starting workers...')

            for key in functypes:
                worker_function, number_workers, joinval = functypes[key]
                work_processes[key] = []
                for i in range(int(number_workers)):
                    t = mp.Process(target=worker_function, args=(i, q, retq, matq, lock))
                    t.daemon = True
                    work_processes[key].append(t)

            for key in work_processes:
                for process in work_processes[key]:
                    logging.info('Starting process')
                    process.start()

            for data in data_generator:
                logging.info('Putting data')
                q.put(data)

            logging.info('Spooling data: OK')

            for key in work_processes:
                for i in range(2):
                    for _ in work_processes[key]:
                        q.put(None)
            logging.info('Joining processes')

            for key in work_processes:
                worker_function, number_workers, joinval = functypes[key]

                if joinval:
                    for process in work_processes[key]:
                        process.join()

            for key in work_processes:
                worker_function, number_workers, joinval = functypes[key]

                if not joinval:
                    for _ in work_processes[key]:
                        retq.put(None)
                        matq.put(None)

            for key in work_processes:
                worker_function, number_workers, joinval = functypes[key]

                if not joinval:
                    for process in work_processes[key]:
                        process.join()

            gc.collect()
            print('DONE!')

        mp_with_timeout(functypes=functype_dict, data_generator=self.generates_dataframes(100))
        return self.hashes_path, self.lshforestpath, self.mat_path


def arg_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('--outpath', help='name of the db', type=str)
    parser.add_argument('--dbtype', help='preconfigured taxonomic ranges', type=str)
    parser.add_argument('--OMA', help='use oma data ', type=str)
    parser.add_argument('--nthreads', help='nthreads for multiprocessing', type=int)
    parser.add_argument('--outfolder', help='folder for storing hash, db and tree objects', type=str)
    parser.add_argument('--verbose', help='print verbose output', type=bool)

    args = parser.parse_args()

    return args


def main(args=None):
    if args is None:
        args = arg_parser()

    dbdict = {
        'all': {'taxfilter': None, 'taxmask': None},
        'plants': {'taxfilter': None, 'taxmask': 33090},
        'archaea': {'taxfilter': None, 'taxmask': 2157},
        'bacteria': {'taxfilter': None, 'taxmask': 2},
        'eukarya': {'taxfilter': None, 'taxmask': 2759},
        'protists': {'taxfilter': [2, 2157, 33090, 4751, 33208], 'taxmask': None},
        'fungi': {'taxfilter': None, 'taxmask': 4751},
        'metazoa': {'taxfilter': None, 'taxmask': 33208},
        'vertebrates': {'taxfilter': None, 'taxmask': 7742},
    }

    if 'outpath' in args:
        dbname = args.outpath
    else:
        raise Exception(' please give your profile an output path with the --outpath argument ')

    if args.dbtype:
        taxfilter = dbdict[args.dbtype]['taxfilter']
        taxmask = dbdict[args.dbtype]['taxmask']
    else:
        taxfilter=None
        taxmask=None

    nperm = 256

    if args.OMA:
        omafile = args.OMA
    else:
        raise Exception(' please specify input data ')

    if args.verbose:
        verbose = args.verbose
    else:
        verbose = False

    threads = 4

    if args.nthreads:
        threads = args.nthreads

    start = time.time()

    with open_file(omafile, mode="r") as h5_oma:
        logging.info('Starting LSH builder')
        lsh_builder = LSHBuilder(h5_oma=h5_oma,
                                 saving_name=dbname,
                                 verbose=verbose,
                                 numperm=nperm,
                                 taxfilter=taxfilter,
                                 taxmask=taxmask
                                 )
    lsh_builder.run_pipeline(threads)
    print(time.time() - start)
    print('DONE')


if __name__ == '__main__':
    args = argparse.Namespace(outpath='out', dbtype='eukarya', OMA='data/OmaServer.h5', verbose=True, nthreads=8)
    main(args)
