

#!/usr/bin/env python
# coding: utf-8
import sys
sys.path.append('../../..')
#sys.path.append( '/home/cactuskid13/miniconda3/pkgs/')
print(sys.path)

from pyprofiler.utils import hashutils
import ete3
import random
from pyprofiler.utils import config_utils
import pyprofiler.utils.goatools_utils as goa
import pyprofiler.utils.hashutils as hashutils
from matplotlib import pyplot as plt
import numpy as np
import pyprofiler.profiler as profiler
import pandas as pd
import itertools
import dask
import warnings
import redis


from sklearn.metrics import precision_recall_curve
import scipy
from dask import dataframe as dd
import pickle


path = '/work/FAC/FBM/DBC/cdessim2/default/dmoi/'

with open(path + 'datasets/STRING/' + 'gold_standard_profiles.pkl' , 'rb' )as profiles_out:
    stringprofiles = pickle.loads(profiles_out.read())
string_df = pd.DataFrame.from_dict(stringprofiles , orient='index')
#make the profiles for this small set of HOGs
for i, key in enumerate(stringprofiles):
    if i < 10:
        print(key,stringprofiles[key])
#now we have profiles for all HUMAP and COG interactions
#String has interactions from each COG in different species.
#We need a way to check for the presence of interaction within a species for a COG
#for this we will create a bloom filter with all the interactions between our cogs
calc_filter = True

if calc_filter == True:
    from dask.distributed import fire_and_forget
    from dask.distributed import Client, Variable , Queue , Lock ,LocalCluster
    from dask_jobqueue import SLURMCluster
    from dask.distributed import  utils_perf
    from dask.distributed import Client, LocalCluster
    import dask
    import redis
    from bloom_filter2 import BloomFilter
    import lzma
    from dask import dataframe as dd
    distributed = True
if calc_filter == True:
    #using distributed computation on a slurm cluster here. This is my particular config. You will need to alter this: https://distributed.dask.org/en/stable/
    if distributed == True:
        NCORE = 4
        print('deploying cluster')
        cluster = SLURMCluster(
            #change theses settings for your cluster
            walltime='4:00:00',
            n_workers = NCORE,
            cores=NCORE,
            processes = NCORE,
            interface='ib0',
            memory="120GB",
            env_extra=[
            
            path + 'miniconda/etc/profile.d/conda.sh',
            'conda activate ML2'
            ],
            #scheduler_options={'interface': 'ens2f0' },
            #if gpu node
            scheduler_options={'interface': 'ens1f1np1' },
            #extra=["--lifetime", "3h55m", "--lifetime-stagger", "4m"]
        )
        print(cluster.job_script())
    else:
        cluster = LocalCluster()
        client = Client(cluster)
if calc_filter == True:
    if distributed == True:
        print(cluster)
        cluster.scale(jobs = 100)
        print(cluster.dashboard_link)
        client = Client(cluster , timeout='450s' , set_as_default=True )
#find which species each of the cogs has an interaction in
if calc_filter == True:
    link_df = dd.read_csv(path + 'datasets/STRING/protein.links.full.v11.5.txt',  blocksize=75e6 , header = 0, sep = ' ')
    print(link_df)

#compute bloom filters for protein pairs
@dask.delayed
def mapcogs(df ):
    #you need a redis server running on your cluster for this to work. change your ip, port and db number accordingly
    rdb = redis.Redis(host='10.202.12.165', port=6379, db=0)
    if type( df ) == tuple:
        df = df[0]
    protlist1 = list(df.protein1.map(lambda x:str(x).strip()))
    protlist2 = list(df.protein2.map(lambda x:str(x).strip()))
    protlist = list(set(protlist1+protlist2))
    data = rdb.mget(protlist)
    mapper = dict(zip(protlist, data) )
    df['COG1'] = df.protein1.map(mapper)
    df['COG2'] = df.protein2.map(mapper)
    df = df.dropna()
    df['COG1'] = df.COG1.map(lambda x:str(x).replace("b",'').replace("'",'').strip() )
    df['COG2'] = df.COG2.map(lambda x:str(x).replace("b",'').replace("'",'').strip() )
    df['species'] = df.protein1.map(lambda x:x.split('.')[0])
    df['coglinks'] = df.COG1 + '_' + df.COG2 + '_' + df.species
    ret = set(df.coglinks.unique())
    return ret

@dask.delayed
def return_filter(coglinks, verbose = True):
    if type( coglinks ) == tuple:
        coglinks = coglinks[0]
    b=BloomFilter(max_elements=10**8, error_rate=0.001 ,start_fresh = True)
    for p in coglinks:
        b.add( p )
    retlen = len(coglinks)
    return   b , retlen

@dask.delayed
def sumfilter(f1,f2, total ):
    if type( f1 ) == tuple:
        f1 = f1[0]
    if type( f2 ) == tuple:
        f2 = f2[0]
    f3 = f1.__ior__(f2)
    return f3 , total

def treesum(totalfilter):
    print(len(totalfilter))
    while len(totalfilter)>1:
        next_round= []
        for i in range(0,len(totalfilter),2):
            if i+1 < len(totalfilter):
                next_round.append( sumfilter( totalfilter[i][0] , totalfilter[i+1][0] , totalfilter[i][1]+totalfilter[i+1][1]  ) )
        if len(totalfilter) % 2 !=0:
            next_round.append(totalfilter[-1])
        totalfilter = next_round
        print(len(totalfilter))
    return totalfilter

if calc_filter == True:
    b=BloomFilter(max_elements=10**8, error_rate=0.001 ,start_fresh = True)
    partitions  = link_df.to_delayed()
    print('map cogs')
    res1 = [ mapcogs(p) for p in partitions ]
    print('done')
    print('make filters')
    res2 = [ return_filter(p) for p in res1 ]
    finals =[]
    for chunk in range(int(len(res2)/1024)+1):
        print(chunk*1024)
        res3 = res2[chunk*1024:(chunk+1)*1024]
        res4 = treesum(res3)
        res4 = dask.compute(res4)
        print(res4)
        finals.append(res4[0])

    with open('bloomfinal_big.pkl' , 'wb' ) as finalout:
        finalout.write(pickle.dumps(finals))

if calc_filter == True:
    with open('bloomfinal_big.pkl' , 'wb' ) as finalout:
        finalout.write(pickle.dumps(finals))
        
with open('bloomfinal_big.pkl' , 'rb' ) as finalout:
    resfinal = pickle.loads(finalout.read()) 
print(resfinal)

