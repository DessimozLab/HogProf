#import keras
import profiler
import os
import sys

from validation import validation_semantic_similarity
from keras.models import Sequential ,model_from_json
from keras.layers import Dense, Activation , Dropout
from keras import optimizers , regularizers , constraints

import pandas as pd
from utils import config_utils , hashutils
import numpy as np
import argparse
import glob
import time as t
import numpy
import pickle
import random
from utils import hashutils
np.random.seed(0)
random.seed(0)
import gc


"""

use this script in conjuction with a labelled protein interaction datasets (string, corum, mips etc).

You can optimize the profile weights at diff taxonomic levels to give a
better prediction of protein interaction with the weighted jaccard score.

A positive and negative dataset need to be generated. This dataset can be from one
organism or several different organisms as long as the truth values are known.


"""


def calculate_x(row):
    mat_x1 = row.mat_x
    mat_x2 = row.mat_y

    ret1 = np.zeros(mat_x1.shape)
    ret2 = np.zeros(mat_x2.shape)
    diff = mat_x1 - mat_x2
    matsum = mat_x1 + mat_x2
    ret1[np.where(diff != 0 ) ] = 1
    ret2[np.where(matsum == 2 ) ] = 1
    return list(np.hstack([ret1,ret2]))


    #generate dataframes w sem sim and profile distances

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument("-hognetcsv", help="csv for training", type =str)
    parser.add_argument("-epochs", help="number of epochs to train", type=int)
    parser.add_argument("-savedir", help="save directory for model", type=str)
    parser.add_argument("-chunksize", help="num hog pairs to analyze at once", type=str)
    parser.add_argument("-overwrite", help="overwrite model", type=str)
    parser.add_argument("-name", help="name of database", type=str)
    parser.add_argument("-traintest", help="train test dataset that has been prebuilt", type=str)
    parser.add_argument("-ntrain", help="n profiles to train on ", type=str)
    parser.add_argument("-ntest", help="n profiles to test on ", type=str)

    args = vars(parser.parse_args(sys.argv[1:]))

    try:
        saving_path  = config_utils.datadir + args['name']

        print(saving_path)

    except:
        raise Exception('please specidfy db name')

    hashes_path = saving_path + 'hashes.h5'
    lshpath = saving_path + 'newlsh.pkl'
    lshforestpath = saving_path + 'newlshforest.pkl'
    mat_path = saving_path+ 'hogmat.h5'

    #load hogs dataset from paper
    #hog names need to be mapped to interactors
    if args['savedir']:
        savedir = args['savedir']
        if not os.path.exists(savedir):
            os.makedirs(savedir)
    else:
        savedir = './'
        # load json and create model

    with open( config_utils.datadir + 'taxaIndex.pkl', 'rb')as taxain:
        taxaIndex = pickle.loads( taxain.read() )
    hogmat_size = 3 *  len(taxaIndex)

    if args['hognetcsv']:
        if args['ntrain']:
            ntrain = int(args['ntrain'])
        else:
            ntrain = 3000
        if args['ntest']:
            ntest = int(args['ntest'])
        else:
            ntest = 1000


        csvs = glob.glob(args['hognetcsv'] )
        print(csvs)
        df = pd.concat ( [ pd.read_csv(csv) for csv in csvs] )
        csvs = glob.glob(args['hognetcsv'] )
        print(csvs)
        df = pd.concat ( [ pd.read_csv(csv) for csv in csvs] )

        #shuffle
        df['HogFamA'] = df.HogA.map(hashutils.hogid2fam)
        df['HogFamB'] = df.HogB.map(hashutils.hogid2fam)
        df = df.sample(frac =1)
        msk = np.random.rand(len(df)) < 0.90
        #split
        print(df)

        traindf = df.iloc[:ntrain,:]
        testdf = df.iloc[ntrain:ntest,:]
        validation = df.iloc[ntest:,:]
        validation.to_csv( savedir + 'validationset.csv')
        print( traindf)
        print( testdf)
        traintest = None
        if 'forest' in args and 'hashes' in args:
            p = profiler.Profiler( hashes = args['hashes'], forest = args['forest'] , oma = True)
        else:
            p = profiler.Profiler( lshforestpath, hashes_path ,  oma = True)

        chunksize = 25
        tstart= t.time()
        gendata= p.retmat_mp(traindf, nworkers = 2, chunksize=50)
        xtotal = []
        ytotal = []
        print('generate data for training')
        for i in range(int(ntrain/ chunksize ) ):
            X,y = next(gendata)
            xtotal.append(X)
            ytotal.append(y)
        xtotalmat = np.vstack(xtotal)
        ytotalmat = np.hstack(ytotal)
        xtotal = []
        ytotal = []
        gendata= p.retmat_mp(testdf, nworkers = 25, chunksize=50)
        print('generate data for testing')
        for i in range(int(ntest/ chunksize ) ):
            X,y  = next(gendata)
            xtotal.append(X)
            ytotal.append(y)
        xtesttotalmat = np.vstack(xtotal)
        ytesttotalmat = np.hstack(ytotal)

        with open( savedir + 'traintest.pkl', 'wb') as traintestout:
            traintestout.write( pickle.dumps( [xtotalmat, ytotalmat, xtesttotalmat, ytesttotalmat ] ) )

    elif args['traintest']:
        #load precomputed profiles to avoid calculating the again
        traintest = args['traintest']
    else:
        raise Exception('please specify input data')

    if args['epochs']:
        eps =  args['epochs']
    else:
        eps = 10
    if args['overwrite'] == 'True':
        overwrite = True
    else:
        overwrite = False



    if args['chunksize']:
        chunksize = args['chunksize']
    else:
        chunksize = 25


    # 3 events, diff and union of matrows
    if overwrite ==False & os.path.isfile(savedir + 'model.json') and  os.path.isfile(savedir + 'model.h5'):
        json_file = open(savedir + 'model.json', 'r')
        loaded_model_json = json_file.read()
        json_file.close()
        model = model_from_json(loaded_model_json)
        # load weights into new model
        model.load_weights(savedir + "model.h5")
        print("Loaded model from disk")
    else:
        print('new model')
        layers = []
        layers.append(Dense( 1, input_dim= hogmat_size  , activation='relu' , use_bias=True))
        layers.append(Dense(1  , activation='sigmoid',  kernel_initializer='random_uniform' )     )

        #layers.append(Dense( 1, activation='softmax' )     )
        #layers.append( Dropout(.5 , noise_shape=None, seed=None))
        model = Sequential(layers)

    #sgd = optimizers.SGD(lr= 1, momentum=0.1, decay=0.01, nesterov=False)
    #sgd = optimizers.Adagrad(lr=0.01, epsilon=.01, decay=0.01)
    sgd = optimizers.Adadelta(lr=.10, rho=0.095, epsilon=None, decay=0.0)
    #sgd = optimizers.RMSprop(lr=1, rho=0.9, epsilon=None, decay=0.0)
    model.compile(loss='binary_crossentropy', optimizer=sgd, metrics=['accuracy'])
    tstart = t.time()

    print('loading dataset')
    with open( traintest , 'rb') as traintestin:
        xtotalmat, ytotalmat, xtesttotalmat, ytesttotalmat  =  pickle.loads( traintestin.read() )
    print( 'done')
    print('training')
    metrics = model.fit(x=xtotalmat  , y=ytotalmat, batch_size= 500 , epochs=2000, verbose=1 )
    print('done')
    print('testing')
    metrics = model.evaluate(x = xtotalmat , y = ytotalmat)
    print(metrics)
    print('done')
    print('saving')
    model_json = model.to_json()
    with open(savedir + "model_nobias.json", "w") as json_file:
        json_file.write(model_json)
    # serialize weights to HDF5
    model.save_weights(savedir + "model_nobias.h5")
    print("Saved model to disk")
    tstart = t.time()
    print('DONE')
