from lshbuilder import LSHBuilder
import profiler
from validation import validation_semantic_similarity
from utils import config_utils
from functools import partial
import pickle
import functools
import itertools
from bayes_opt import BayesianOptimization
from utils import hashutils
import numpy as np
import random
import multiprocessing as mp
import time
import gc

from tables import *
import argparse, sys
import os

from bayes_opt.observer import JSONLogger
from bayes_opt.event import Events
from bayes_opt.util import load_logs
import gc

import glob
#set rand seed
np.random.seed(0)
random.seed(0)


parser=argparse.ArgumentParser()
parser.add_argument('--lw', help='loss weight',type = float)
parser.add_argument('--pw', help='presence weight',type = float)
parser.add_argument('--dw', help='duplication weight' , type = float)

parser.add_argument('--ll', help='loss lambda', type = float)
parser.add_argument('--pl', help='presence lambda' , type = float)
parser.add_argument('--dl', help='duplication lambda', type = float )

parser.add_argument('--db', help='DB type' , type=str)
parser.add_argument('--dir', help='save dir')

args=parser.parse_args()


dbdict = {
'plants': { 'taxfilter': None , 'taxmask': 33090 },
'all': { 'taxfilter': None , 'taxmask': None },
'archaea':{ 'taxfilter': None , 'taxmask': 2157 },
'bacteria':{ 'taxfilter': None , 'taxmask': 2 },
'eukarya':{ 'taxfilter': None , 'taxmask': 2759 },
'protists':{ 'taxfilter': [2 , 2157 , 33090 , 4751, 33208] , 'taxmask':None },
'fungi':{ 'taxfilter': None , 'taxmask': 4751 },
'metazoa':{ 'taxfilter': None , 'taxmask': 33208 },
}

def profiling_error( db , taxfilter, tax_mask, lossweight , presenceweight, dupweight, loss_lambda , presence_lambda , dupl_lamba,  hoglist , val = None, compile = True , dir = None):
    print('compiling' + db)
    #record param settings
    #compile lsh
    parastr = 'lw'+str(lossweight)+ 'pw'+str(presenceweight)+ 'dw'+ str(dupweight)+ 'll'+str(loss_lambda)+ 'pl'+ str(presenceweight) +'dl' + str(dupl_lamba)
    #print(parastr)
    startdict={'presence':presenceweight, 'loss':lossweight, 'dup':dupweight}
    lambdadict={'presence':presence_lambda, 'loss':loss_lambda, 'dup':dupl_lamba}
    if compile == True:
        with open_file(config_utils.omadir + 'OmaServer.h5', mode="r") as h5_oma:
            lsh_builder = LSHBuilder(h5_oma, saving_folder= dir , saving_name=db, numperm = 512,
            treeweights= None , taxfilter = taxfilter, taxmask= tax_mask , lambdadict= lambdadict, start= startdict)
            hashes, forest , mat = lsh_builder.run_pipeline()
            #hashes, forest, lshpath =lsh_builder.run_pipeline()
        print( 'done compiling')
    else:
        hashes = dir + 'hashes.h5'
        forest = dir + 'newlshforest.pkl'
        mat = dir + 'hogmat.h5'
    print('query DB and calculate error')
    print('load profiler')
    p = profiler.Profiler(lshforestpath = forest, hashes_h5=hashes, mat_path= None , nsamples = 512)

    print('done')
    print('loading validation')
    if val is None:
        if not os.path.isfile(config.datadir + 'val.pkl'):
            folder = config_utils.datadir + 'GOData/'
            val = validation_semantic_similarity.Validation_semantic_similarity( folder + 'go-basic.obo' ,
                folder + 'goframe.pkl' , folder + 'oma-go.txt' , config_utils.omadir + 'OmaServer.h5' , folder + 'termcounts.pkl' )
            with open(config.datadir + 'val.pkl' , 'wb')as valout:
                valout.write(pickle.dumps(val))
        else:
            with open(config.datadir + 'val.pkl' , 'rb')as valout:
                val = pickle.loads(valout.read())
    print( 'done')
    print('testing db')
    if not hoglist:
        #sample random hogs
        hoglist = list(np.random.randint(0, high=610000, size=200, dtype='l'))
        hoglist = [ hashutils.fam2hogid(s)  for s in hoglist]
    scores = {}
    retq = mp.Queue()
    lock = mp.Lock()
    timelimit = 100
    for i,hog in enumerate(hoglist):
        res = p.hog_query( hog_id = hog , k = 20)
        res = set([ hashutils.fam2hogid(r) for r in res]+[hog])
        #write loop for sem sim check with timeout here
        processes = {}
        for combo in itertools.combinations(res,2):
            processes[combo] = {'time':time.time() , 'process': mp.Process( target = val.semantic_similarity_score_mp , args = (combo[0],combo[1],retq , lock)  ) }
            processes[combo]['process'].start()
            while len(processes)> mp.cpu_count()/4:
                time.sleep(.01)
                for c in processes:
                    if processes[c]['time']>timelimit or processes[c]['process'].exitcode is not None:
                        #print( c[0] +':' + c[1] + ' done')
                        processes[c]['process'].terminate()
                        gc.collect()
                        del(processes[c])
                        break
        while len(processes)> 0:
            time.sleep(.01)
            for c in processes:
                if processes[c]['time']>timelimit or processes[c]['process'].exitcode is not None:
                    processes[c]['process'].terminate()
                    if rocesses[c]['time']>timelimit:
                        print('timeout')
                    gc.collect()
                    del(processes[c])
                    break
        hogsemsim = {}

        while retq.empty() == False:
            combo,semsim = retq.get()
            print(combo)
            print(semsim)
            hogsemsim[combo]=semsim

        scores.update( { combo: {'query_num':i, 'hog_sem_sim_normalize': hogsemsim[combo][0] , 'hog_sem_sim': hogsemsim[combo][1]
        , 'hog_jaccard_sim' : p.hog_v_hog(combo[0],combo[1])
        }  for combo in itertools.combinations(res,2)  if combo in hogsemsim } )
        print(scores)
    resdf = pd.DataFrame.from_dict( scores, orient = 'index')
    resdf.to_csv( dir + 'resdf' + '.csv')
    #take positive information values
    semsim_mean = resdf[resdf.hog_sem_sim >0].hog_sem_sim_normalize.mean()
    print(semsim_mean)
    print('done')
    return semsim_mean


if __name__ == '__main__':
    print(sys.argv)
    args = vars(parser.parse_args(sys.argv[1:]))

    db = args['db']
    dir = args['dir']
    if not os.path.exists(dir):
        os.makedirs(dir)
    error = functools.partial( profiling_error , db=db , taxfilter = dbdict[db]['taxfilter'], tax_mask = dbdict[db]['taxmask'],  hoglist =None , dir = args['dir'])
    #get error for the first point with all weights at 1
    bo = BayesianOptimization( f = error ,  pbounds = {'lossweight': (0, 1),
                                                    'presenceweight': (0, 1),
                                                    'dupweight':(0,1),
                                                    'loss_lambda':(-1,1),
                                                    'presence_lambda':(-1,1),
                                                    'dupl_lamba':(-1,1)
                                                    } ,
                                verbose = 2,
                                random_state = 0,
                                                )

    if len(glob.glob("./logs"+db+".json") )>0:
        load_logs(bo, logs=["./logs"+db+".json"])
    logger = JSONLogger(path="./logs"+db+".json")
    bo.subscribe(Events.OPTMIZATION_STEP, logger)
    try:
        if args['lw'] and args['pw'] and args['dw']:
            #try specific points
            lw = args['lw']
            pw = args['pw']
            dw = args['dw']
            ll = 0
            pl = 0
            dl = 0

            try:
                if args['ll'] and args['pl'] and args['dl']:
                    ll = args['ll']
                    pl = args['pl']
                    dl = args['dl']
            except:
                pass

            bo.probe(
            params={'lossweight':  lw,
                    'presenceweight': pw,
                    'dupweight': dw,
                    'loss_lambda':ll,
                    'presence_lambda':pl,
                    'dupl_lamba':dl
                    },
            lazy=True,
            )
    except:
        pass

    #else use loaded points to try new interesting ones
    bo.maximize(init_points=2, n_iter=15, kappa=2)
    #save the friggin result
    with open( savedir+'bayesopt.pkl', mode='wb', buffering=None) as bayesout:
        bayesout.write(  pickle.dumps(bo, -1))
    print('DONE')
    print(bo.res['max'])
