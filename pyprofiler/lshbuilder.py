from tables import *
import functools
import multiprocessing as mp
import pandas as pd
import time as t
import pickle
from datasketch import MinHashLSH,   MinHashLSHForest , WeightedMinHashGenerator
from datetime import datetime
import h5py
import json
import ete3
import gc
from utils import  config_utils, pyhamutils, hashutils , files_utils
import numpy as np
import random

random.seed(0)
np.random.seed(0)

class LSHBuilder:
    """


    This class contains the stuff you need to make a phylogenetic profiling database with input orthxml files and a taxonomic tree

    You must either input an OMA hdf5 file or an ensembl tarfile containing orthoxml file with orthologous groups.

    You can provide a species tree or use the ncbi taxonomy with a list of taxonomic codes for all the species in your db

    """
    def __init__(self, tarfile_ortho = None,  h5_oma = None, taxa = None, masterTree = None, saving_name=None ,   numperm = 256,  treeweights= None , taxfilter = None, taxmask= None ,  verbose = False):
        if h5h5_oma:
            from pyoma.browser import db
            self.h5OMA = h5_oma
            self.db_obj = db.Database(h5_oma)
            self.oma_id_obj = db.OmaIdMapper(self.db_obj)
        elif tarfile_ortho:
            self.tar = tarfile_ortho
            self.h5OMA = None
            self.db_obj = None
            self.oma_id_obj = None
        else:
            raise Exception( 'please specify input data. Either a tarfile of orthxml files or OMA hdf5 ' )
        self.tax_filter = taxfilter
        self.tax_mask = taxmask
        self.verbose = verbose
        self.saving_folder = saving_folder
        self.datetime = datetime
        self.date_string = "{:%B_%d_%Y_%H_%M}".format(datetime.now())
        self.saving_name= saving_name
        if saving_name:
            self.saving_path =config_utils.datadir + saving_name +'/'
            if not os.path.isdir(self.saving_path):
                os.mkdir(path=self.saving_path)
        else:
            self.saving_path = config_utils.datadir + self.date_string +'/'
            if not os.path.isdir(self.saving_path):
                os.mkdir(path=self.saving_path)

        if masterTree is None:
            if h5h5_oma:
                genomes = pd.DataFrame(h5_oma.root.Genome.read())["NCBITaxonId"].tolist()
                genomes = [ str(g) for g in genomes]
                taxa = genomes + [ 131567, 2759, 2157, 45596 ]+[ taxrel[0] for taxrel in  list(h5_oma.root.Taxonomy[:]) ]  + [  taxrel[1] for taxrel in list(h5_oma.root.Taxonomy[:]) ]
                self.tree_string , self.tree_ete3 = files_utils.get_tree(taxa=taxa , savename =saving_name )
            elif taxa:
                with open(taxa, 'r') as taxin:
                    taxlist = [ int(line) for line in taxin ]
                self.tree_string , self.tree_ete3 = files_utils.get_tree(taxa=taxlist , savename =saving_name )
            else:
                raise Exception( 'please specify either a list of taxa or a tree' )
        elif mastertree:
            with open( masterTree , 'wb') as pklin:
                self.tree_ete3 = pickle.loads(pklin.read())
                self.tree_string = self.tree_ete3.write(format=1)

        self.taxaIndex, self.reverse = files_utils.generate_taxa_index(self.tree_ete3 , self.tax_filter, self.tax_mask)

        with open( config_utils.saving_path + 'taxaIndex.pkl', 'wb') as taxout:
            taxout.write( pickle.dumps(self.taxaIndex))

        self.numperm = numperm
        if treeweights is None:
            #generate all ones
            self.treeweights = hashutils.generate_treeweights(self.tree_ete3  , self.taxaIndex , taxfilter, taxmask)
        else:
            self.treeweights = treeweights

        wmg = WeightedMinHashGenerator(3*len(self.taxaIndex), sample_size=numperm, seed=1)
        with open( self.saving_path +saving_name + 'wmg.pkl', 'wb') as taxout:
            taxout.write( pickle.dumps(self.taxaIndex))

        self.wmg = wmg
        self.HAM_PIPELINE = functools.partial(pyhamutils.get_ham_treemap_from_row, tree=self.tree_string )
        self.HASH_PIPELINE = functools.partial(hashutils.row2hash , taxaIndex=self.taxaIndex  , treeweights=self.treeweights , wmg=wmg )
        if self.h5OMA:
            self.READ_ORTHO = functools.partial(pyhamutils.get_orthoxml_oma, db_obj=self.db_obj)
        elif self.tar:
            self.READ_ORTHO = pyhamutils.get_orthoxml
        self.hashes_path = self.saving_path + 'hashes.h5'
        self.lshpath = self.saving_path + 'newlsh.pkl'
        self.lshforestpath = self.saving_path + 'newlshforest.pkl'
        self.mat_path = self.saving_path+ 'hogmat.h5'
        self.columns = len(self.taxaIndex)

    def load_one(self, fam):
        #test function to try out the pipeline on one orthoxml
        ortho_fam = self.READ_ORTHO(fam)
        pyham_tree = self.HAM_PIPELINE([fam, ortho_fam])
        hog_matrix,weighted_hash = hashutils.hash_tree(pyham_tree , self.taxaIndex , self.treeweights , self.wmg)
        return ortho_fam , pyham_tree, weighted_hash,hog_matrix

    def generates_dataframes(self, size=100, minhog_size=3, maxhog_size=None ):
        families = {}
        start = -1
        if self.h5OMA:
            self.groups  = self.h5OMA.root.OrthoXML.Index
            self.rows = len(self.groups)
            for i, row in enumerate(self.groups):
                if i > start:
                    fam = row[0]
                    ortho_fam = self.READ_ORTHO(fam)
                    hog_size = ortho_fam.count('<species name=')
                    if (maxhog_size is None or hog_size < maxhog_size) and (minhog_size is None or hog_size > minhog_size):
                        families[fam] = {'ortho': ortho_fam}
                    if len(families) > size:
                        pd_dataframe = pd.DataFrame.from_dict(families, orient='index')
                        pd_dataframe['Fam'] = pd_dataframe.index
                        yield pd_dataframe
                        families = {}
            pd_dataframe = pd.DataFrame.from_dict(families, orient='index')
            pd_dataframe['Fam'] = pd_dataframe.index
            yield pd_dataframe
            print('last dataframe sent')
            families = {}


        elif self.tar:
            groupfiles=glob.glob(self.tar)
            groups = []
            for tarfile in groupfiles:
                with tarfile.open(tarfile, "r:gz") as tar:
                    for member in tar.getmembers():
                        f=tar.extractfile(member)
                        oxml = ET.parse(f)
                        for input in oxml.iter():
                            ortho_fam = ET.tostring( next(oxml.iter()), encoding='utf8', method='xml' ).decode() )
                            hog_size = ortho_fam.count('<species name=')
                            if (maxhog_size is None or hog_size < maxhog_size) and (minhog_size is None or hog_size > minhog_size):
                                families[member] = {'ortho': ortho_fam}
                    tar.close()


    def universe_saver(self, i, q, retq, matq,univerq, l):
        #only useful to save all prots within a taxonomic range as db is being compiled
        allowed = set( [ n.name for n in self.tree_ete3.get_leaves() ] )
        with open(self.saving_path+'universe.txt') as universeout:
            while True:
                prots = univerq.get()
                for row in df.iterrows():
                    for ID in row.prots.tolist():
                        universeout.write(ID)
                else:
                    print('Universe saver done' + str(i))
                    break

    def worker(self, i, q, retq, matq, l):
        if self.verbose == True:
            print('worker init ' + str(i))
        while True:
            df = q.get()
            if df is not None :
                df['tree'] = df[['Fam', 'ortho']].apply(self.HAM_PIPELINE, axis=1)
                df[['hash','rows']] = df[['Fam', 'tree']].apply(self.HASH_PIPELINE, axis=1)
                retq.put(df[['Fam', 'hash']])
                #matq.put(df[['Fam', 'rows']])
            else:
                if self.verbose == True:
                    print('Worker done' + str(i))
                break

    def saver(self, i, q, retq, matq, l):
        print_start = t.time()
        save_start = t.time()
        global_time = t.time()
        chunk_size = 100
        count = 0
        forest = MinHashLSHForest(num_perm=self.numperm)

        taxstr = ''
        if self.tax_filter is None:
            taxstr = 'NoFilter'
        if self.tax_mask is None:
            taxstr+= 'NoMask'
        else:
            taxstr = str(self.tax_filter)
        dataset_name = self.saving_name+'_'+taxstr
        self.errorfile = self.saving_path + 'errors.txt'
        with open(self.errorfile, 'w') as hashes_error_files:
            with h5py.File(self.hashes_path, 'w', libver='latest') as h5hashes:
                datasets = {}
                if dataset_name not in h5hashes.keys():
                    if self.verbose == True:
                        print('creating dataset')
                        print(dataset_name)
                        print('filtered at taxonomic level: '+taxstr)
                    h5hashes.create_dataset(dataset_name+'_'+taxstr, (chunk_size, 0), maxshape=(None, None), dtype='int32')
                    datasets[dataset_name] = h5hashes[dataset_name+'_'+taxstr]
                    if self.verbose == True:
                        print(datasets)
                    h5flush = h5hashes.flush
                print('saver init ' + str(i))
                while True:
                    this_dataframe = retq.get()
                    if this_dataframe is not None:
                        if not this_dataframe.empty:
                            hashes = this_dataframe['hash'].to_dict()
                            print(str(this_dataframe.Fam.max())+ 'fam num')
                            print(str(count) + ' done')
                            hashes = {fam:hashes[fam] if hashes[fam] is not None else print(fam) for fam in hashes }
                            [ forest.add(str(fam),hashes[fam]) for fam in hashes]
                            for fam in hashes:
                                if len(datasets[dataset_name]) < fam + 10:
                                    datasets[dataset_name].resize((fam + chunk_size, len(hashes[fam].hashvalues.ravel())))
                                datasets[dataset_name][fam, :] = hashes[fam].hashvalues.ravel()
                                count += 1
                            if t.time() - save_start > 200:
                                print( t.time() - global_time )
                                forest.index()
                                print(forest.query( hashes[fam] , k = 10 ) )
                                h5flush()
                                save_start = t.time()
                                with open(self.lshforestpath , 'wb') as forestout:
                                    forestout.write(pickle.dumps(forest, -1))
                                if self.verbose == True:
                                    print('save done at' + str(t.time() - global_time))
                        else:
                            print(this_dataframe)
                    else:
                        if self.verbose == True:
                            print('wrap it up')
                        with open(self.lshforestpath , 'wb') as forestout:
                            forestout.write(pickle.dumps(forest, -1))
                        h5flush()
                        if self.verbose == True:
                            print('DONE SAVER' + str(i))
                        break


    def matrix_updater(self, iprocess , q, retq, matq, l):
        save_start = t.time()
        chunk_size = 100
        print('hogmat saver init ' + str(iprocess))
        h5mat = None
        times1 = []
        times2 = []
        frames = []
        with h5py.File(self.mat_path + str(iprocess) + 'h5' , 'w', libver='latest' ) as h5hashes:
            i =0
            while True:
                rows = matq.get()
                if rows is not None:
                    rows = rows.dropna()
                    maxfam = rows.Fam.max()
                    if h5mat is None:
                        h5hashes.create_dataset( 'matrows',(10,block.shape[1]), maxshape=(None, block.shape[1]),chunks=(1, block.shape[1]), dtype='i8')
                        h5mat = h5hashes['matrows']
                    if h5mat.shape[0] < maxfam:
                        h5mat.resize((maxfam+1,block.shape[1]))
                    i+=1
                    frames.append(rows)
                    assign = t.time()
                    index = np.asarray(rows.Fam)
                    block = np.vstack(rows.rows)
                    h5mat[index,:]= block

                    times1.append(t.time()-assign)
                    if len(times1)>10:
                        times1.pop(0)
                        print(np.mean(times1))
                    h5hashes.flush()
                else:
                    h5hashes.flush()
                    break
        print('DONE MAT UPDATER' + str(i))

    def run_pipeline(self):
        functype_dict = {'worker': (self.worker, int(1*mp.cpu_count()/2), True), 'updater': (self.saver, 1, False),
                         'matrix_updater': (self.matrix_updater, 0, False) }
        self.mp_with_timeout(functypes=functype_dict, data_generator=self.generates_dataframes(100))
        return self.hashes_path, self.lshforestpath , self.mat_path

    @staticmethod
    def mp_with_timeout(functypes, data_generator):
        work_processes = {}
        update_processes = {}
        lock = mp.Lock()
        cores = mp.cpu_count()
        q = mp.Queue(maxsize=cores * 10)
        retq = mp.Queue(maxsize=cores * 10)
        matq = mp.Queue(maxsize=cores * 10)
        work_processes = {}
        print('start workers')
        for key in functypes:
            worker_function, number_workers, joinval = functypes[key]
            work_processes[key] = []
            for i in range(int(number_workers)):
                t = mp.Process(target=worker_function, args=(i, q, retq, matq, lock))
                t.daemon = True
                work_processes[key].append(t)
        for key in work_processes:
            for process in work_processes[key]:
                process.start()
        count = 0
        for data in data_generator:
            q.put(data)
        print('done spooling data')
        for key in work_processes:
            for i in range(2):
                for _ in work_processes[key]:
                    q.put(None)
        print('joining processes')
        for key in work_processes:
            worker_function, number_workers , joinval = functypes[key]
            if joinval == True:
                for process in work_processes[key]:
                    process.join()

        for key in work_processes:
            worker_function, number_workers, joinval = functypes[key]
            if joinval == False:
                for _ in work_processes[key]:
                    retq.put(None)
                    matq.put(None)

        for key in work_processes:
            worker_function, number_workers , joinval = functypes[key]
            if joinval == False:
                for process in work_processes[key]:
                    process.join()
        gc.collect()
        print('DONE!')



if __name__ == '__main__':

    parser=argparse.ArgumentParser()
    parser.add_argument('--taxweights', help='load optimised weights from keras model',type = str)
    parser.add_argument('--taxmask', help='consider only one branch',type = str)
    parser.add_argument('--taxfilter', help='remove these taxa' , type = str)
    parser.add_argument('--name', help='name of the db', type = str)
    parser.add_argument('--dbtype', help='preconfigured taxonomic ranges' , type = str)
    parser.add_argument('--OMA', help='use oma data ' , type = str)
    parser.add_argument('--tarfile', help='use tarfile with orthoxml data ' , type = str)
    parser.add_argument('--nperm', help='number of hash functions to use when constructing profiles' , type = str)
    parser.add_argument('--masterTree', help='number of hash functions to use when constructing profiles' , type = str)

    parser.add_argument('--masterTree', help='number of hash functions to use when constructing profiles' , type = str)

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

    args = vars(parser.parse_args(sys.argv[1:]))

    if 'name' in args:
        dbname = args['name']
    else:
        raise Exception(' please give your profile db a name ')
    if 'dbtype' in args:
        taxfilter = dbdict[args['dbtype']]['taxfilter']
        taxmask = dbdict[args['dbtype']]['taxmask']
    if 'taxmask' in args:
        taxfilter = json.loads ( args['taxmask'] )
    if 'taxfilter' in args:
        taxmask = args['taxfilter']

    if 'nperm' in args:
        nperm = args['nperm']
    else:
        nperm = 256

    if 'omafile' in args:
        omafile = args['omafile']
        tarfile = None
    elif 'tarfile' in args:
        tarfile = args['tarfile']
        omafile = None
    elif config_utils.omadir:
        omafile = config_utils.omadir + 'OmaServer.h5'
        tarfile = None
    elif config_utils.tarfile:
        tarfile = config_utils.tarfile
        omafile = None
    else:
        raise Exception(' please specify input data ')

    if 'taxweights' in args:
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

    print('compiling' + dbname)
    if 'masterTree' in args:
        mastertree = args['masterTree']
    else:
        mastertree=None


    import resource

    start = time.time()
    if omafile:

        with open_file( omafile , mode="r") as h5_oma:
            lsh_builder = LSHBuilder(h5_oma = h5_oma,  saving_name=dbname, numperm = nperm ,
            treeweights= weights , taxfilter = taxfilter, taxmask=taxmask , masterTree =mastertree )
            lsh_builder.run_pipeline()
    else:
        lsh_builder = LSHBuilder( tarfile_ortho = tarfile ,  saving_name=dbname, numperm = perm ,
        treeweights= weights , taxfilter = taxfilter, taxmask=taxmask , masterTree =mastertree )
        lsh_builder.run_pipeline()

    r = resource.getrusage(resource.RUSAGE_BOTH).ru_maxrss
    print('maxmem')
    print(r)
    print('time')
    print(time.time() - start)

    print('DONE')
