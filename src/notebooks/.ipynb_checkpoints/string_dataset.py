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


#load profiler
p = profiler.Profiler(lshforestpath = path + 'datasets/OMA/sep2022/all/newlshforest.pkl' , hashes_h5=path +'datasets/OMA/sep2022/all/hashes.h5' , mat_path= None, oma = path + 'datasets/OMA/sep2022/OmaServer.h5', tar= None , nsamples = 256 , mastertree = path + 'datasets/OMA/sep2022/all_master_tree.pkl')
def grabHog(ID, verbose = True):
    try:
        entry = p.db_obj.entry_by_entry_nr(p.db_obj.id_resolver.resolve(ID))
        return entry[4].decode() , entry
    except:
        return np.nan,np.nan


#################begin building the string dataset ###################################

filter_coglinks = True
#filtering down the COGlinks using cutoffs for global score and text mining
if filter_coglinks == True:
    from collections import Counter
    coglink_df = dd.read_csv(path +'datasets/STRING/COG.links.detailed.v11.5.txt', blocksize=25e6 , header = 0, sep = ' ')
    print(coglink_df)
    print(coglink_df.columns)
    dropcols = [ 'cooccurence', 'combined_score' ]
    coglink_df = coglink_df.drop(columns = dropcols)
    coglink_df['score'] = coglink_df.coexpression + coglink_df.experimental +coglink_df.database+ coglink_df.textmining
    #these cutoffs ere found below using jaccard and AUC
    coglink_df= coglink_df[coglink_df.score>1000]
    coglink_df= coglink_df[coglink_df.textmining>500]
    coglink_df = coglink_df.compute()
    coglink_count = Counter(list(coglink_df.group1)+list(coglink_df.group2))
    coglink_df['count1']= coglink_df.group1.map(coglink_count)
    coglink_df['count2']= coglink_df.group2.map(coglink_count)
    #filter input set
    #coglink_df = coglink_df[coglink_df.count1 > 50 ]
    #coglink_df = coglink_df[coglink_df.count2 > 50 ]
    
    print(coglink_df.head() , len(coglink_df))

#map the interacting cogs to the proteins
compute_grabcogs = False
if compute_grabcogs == True:
    grabcogs = set( list(coglink_df.group1.unique()) + list(coglink_df.group2.unique()) )
    grabcogs= list(grabcogs)
    COGmapings_df = dd.read_csv(path +'datasets/STRING/COG.mappings.v11.5.txt', blocksize=25e6 , header = 0, sep = '\t')
    COGmapings_df = COGmapings_df.set_index('orthologous_group')
    COGmapings_df.astype(str)
    COGmapings_df['##protein'].map( lambda x : x.strip() )
    COGmapings_df['species'] = COGmapings_df['##protein'].map( lambda x : x.split('.')[0] )
    COGmapings_df['COG'] = COGmapings_df.index
    COGmapings_df = COGmapings_df.loc[grabcogs]
    COGmapings_df = COGmapings_df.compute()
    print(COGmapings_df.head())
#only take the proteins in our cogs of interest
if compute_grabcogs == True:
    grabprots =list(COGmapings_df['##protein'].unique())
    print(len(grabprots))
    with open(path + 'datasets/STRING/COG.links.detailed.v11.5.txt' + '.grabcogs.txt', 'w') as protsout:
        protsout.write(''.join([ p + '\n' for p in grabcogs ]) )
    with open(path + 'datasets/STRING/COG.mappings.v11.5.txt' + '.grabprots.txt' , 'w') as protsout:
        protsout.write(''.join([ p + '\n' for p in grabprots ]) )
else:
    with open(path + 'datasets/STRING/COG.links.detailed.v11.5.txt' + '.grabcogs.txt', 'r') as protsout:
        grabcogs = [ cog for cog in protsout.readlines()]
    with open(path +'datasets/STRING/COG.mappings.v11.5.txt' + '.grabprots.txt' , 'r') as protsout:
        grabprots = [ prot for prot in protsout.readlines()]
calc_mappers = False
#use redis to store the mapping of proteins to their cogs
if calc_mappers == True:
    rdb = redis.Redis(host='10.202.12.165', port=6379, db=0)
    count = 0
    for i,r in COGmapings_df.iterrows():
        rdb.set(r['##protein'], i)
        count+=1
        if count < 10:
            print(i+'\n',r)
        if count%1000000==0:
            print(count/len(COGmapings_df))
maphogs = False
if maphogs == True:
    #you need to change this for your own rdb configuration
    #this is used later by the distributed computation
    #mapping each string cog to an oma hog by selecting a member of the cog
    rdb = redis.Redis(host='10.202.12.165', port=6379, db=0)
    hogmap = {}
    for i,prot in enumerate(grabprots):
        if i % 100000 == 0 :
            print(i/len(grabprots))
        cog = rdb.get(prot)
        if cog not in hogmap:
            mapped =  grabHog(prot)
            #retry until something maps
            if mapped[0] != np.nan and type(mapped[0]) == str :
                if len(mapped[0])>1 :
                    hogmap[cog] = mapped
    with open('stringhogmap.pkl' , 'wb')as hogmapout:
        hogmapout.write(pickle.dumps(hogmap))
else:
    with open('stringhogmap.pkl' , 'rb')as hogmapout:
        hogmap = pickle.loads(hogmapout.read())
print(len(hogmap))
for i, key in enumerate(hogmap):
    if i < 10:
        print(key, hogmap[key])
#add the HOGs to the COGdf
#grab the corresponding profiles
compile_final_cogdf = True
if compile_final_cogdf == True:
    print(len(coglink_df))
    try:
        coglink_df.group1  = coglink_df.group1.map( lambda x : x.encode())
        coglink_df.group2  = coglink_df.group2.map( lambda x : x.encode())
    except:
        pass
    coglink_df['hog1'] = coglink_df.group1.map(hogmap)
    coglink_df['hog2'] = coglink_df.group2.map(hogmap)
    coglink_df=coglink_df.dropna()
    print(len(coglink_df))
    print(coglink_df.head())
    coglink_df['hogid_1'] = coglink_df['hog1'].map(lambda x:x[0])
    coglink_df['hogid_2'] = coglink_df['hog2'].map(lambda x:x[0])
    coglink_df['fam1'] = coglink_df['hog1'].map( lambda x :   p.hogid2fam(x[1]) )
    coglink_df['fam2'] = coglink_df['hog2'].map( lambda x :   p.hogid2fam(x[1]) ) 
    coglink_df.fam1 = coglink_df.fam1.map(int)
    coglink_df.fam2 = coglink_df.fam2.map(int)
    stringHOGs = set(coglink_df.fam1.unique()).union(set(coglink_df.fam2.unique()))
    print(len(stringHOGs))
    print(coglink_df)
    coglink_df.to_csv('STRINGCOGS2OMAHOGS.csv')
else:
    coglink_df = pd.read_csv('STRINGCOGS2OMAHOGS.csv')
stringPairs = coglink_df
#derive explicit profiles for our hogs of interest in string
calc_hogs_string = True
stringprofiles = {}
if calc_hogs_string == True:
    print('profiles to calclulate',len(stringHOGs))
    for i,fam in enumerate(stringHOGs):
        if i % 100 ==0:
            print(i)
        try:
            prof = p.return_profile_OTF(fam)
            stringprofiles.update(prof)
        except:
            print('err',fam)
if calc_hogs_string == True:
    with open(path + 'datasets/STRING/' + 'gold_standard_profiles.pkl' , 'wb') as profiles_out:
        profiles_out.write(pickle.dumps(stringprofiles))


