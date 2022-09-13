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
from sklearn.metrics import precision_recall_curve
import scipy
from dask import dataframe as dd
import pickle

path = '/work/FAC/FBM/DBC/cdessim2/default/dmoi/'
#lets load a compiled db containing the OMA root HOGs into a profiler oject 
p = profiler.Profiler(lshforestpath = path + 'datasets/OMA/sep2022/all/newlshforest.pkl' , hashes_h5=path +'datasets/OMA/sep2022/all/hashes.h5' , mat_path= None, oma = path + 'datasets/OMA/sep2022/OmaServer.h5', tar= None , nsamples = 256 , mastertree = path + 'datasets/OMA/sep2022/all_master_tree.pkl')
#map ids to OMA HOGs
def grabHog(ID, verbose = True):
    try:
        entry = p.db_obj.entry_by_entry_nr(p.db_obj.id_resolver.resolve(ID))
        return entry[4].decode() , entry
    except:
        return np.nan,np.nan
#filtering down the humap dataset
humap = path + 'datasets/humap_PPI/humap2_ppis_ACC_20200821.pairsWprob'
calc_humap = True
if calc_humap == True:
    #load humap data
    humap_df = pd.read_table(humap, header = None)
    print(humap_df)
    print(len(humap_df))
    humap_df = humap_df[humap_df[2] > 0.03 ]
    print(len(humap_df))
    print('map all hogids')
    mapper = set( list(humap_df[1]) + list(humap_df[0]) )
    mapper = { protid: grabHog(protid) for protid in mapper }
    print('done compiling mapper')
    humap_df['hog1'] = humap_df[1].map(mapper)
    humap_df['hog2'] = humap_df[0].map(mapper)
    humap_df['hogid_1'] = humap_df['hog1'].map(lambda x:x[0])
    humap_df['hogid_2'] = humap_df['hog2'].map(lambda x:x[0])
    humap_df = humap_df.dropna()
    humap_df['fam1'] = humap_df['hog1'].map( lambda x :   p.hogid2fam(x[1]) )
    humap_df['fam2'] = humap_df['hog2'].map( lambda x :   p.hogid2fam(x[1]) ) 
    humap_df = humap_df.dropna()
    humap_df.fam1 = humap_df.fam1.map(int)
    humap_df.fam2 = humap_df.fam2.map(int)
    print(len(humap_df))
    humap_df.to_csv(humap+'hogmapped.csv')
else:
    humap_df = pd.read_csv(humap+'hogmapped.csv')
print(len(humap_df))
humap_pairs = humap_df
#calculating the profile vectors for each HOG in the Humap dataset
calc_hogs_humap = True
if calc_hogs_humap == True:
    allhogs = set([])
    allhogs = allhogs.union( set(humap_df.fam1.unique() ) )
    allhogs = allhogs.union( set(humap_df.fam2.unique() ) )
    allhogs = list(allhogs)
    print(len(allhogs))
if calc_hogs_humap == True:
    profiles = p.retmat_mp_profiles(allhogs)
print(profiles)
save_Hogs_humap = True
if save_Hogs_humap == True:
    with open(humap + 'gold_standard_profiles.pkl' , 'wb') as profiles_out:
        profiles_out.write(pickle.dumps(profiles))