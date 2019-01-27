#import keras
import profiler
import os
import sys

from validation import validation_semantic_similarity

from keras.models import Sequential
from keras.layers import Dense, Activation
import pandas as pd
from utils import config_utils , hashutils
import numpy as np
# MLP for Pima Indians Dataset Serialize to JSON and HDF5
from keras.models import Sequential
from keras.layers import Dense
from keras.models import model_from_json

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
###return profiler and validation obj


def load_valobjs( db , val = False , hashes=None , forest=None ):
    print('compiling' + db)
    p = profiler.Profiler(lshforestpath = forest, hashes_h5=hashes, mat_path= None , nsamples = 256, oma = True )
    print('done')
    print('loading validation')
    if val:
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
    return p, val


def calculate_x(row):
    mat_x1 = row.mat_x
    mat_x2 = row.mat_y

    ret1 = np.zeros(mat_x1.shape)
    ret2 = np.zeros(mat_x2.shape)
    diff = mat_x1 - mat_x2
    matsum = mat_x1 + mat_x2
    ret1[np.where(diff != 0 ) ] = -1
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
    parser.add_argument("-forest", help="lsh forest", type=str)
    parser.add_argument("-hashes", help="hashvals for forest", type=str)
    print(sys.argv)
    args = vars(parser.parse_args(sys.argv[1:]))
    #load hogs dataset from paper
    if args['hognetcsv']:
        csvs = glob.glob(args['hognetcsv'] )
        print(csvs)
        df = pd.concat ( [ pd.read_csv(csv) for csv in csvs] )

    if args['epochs']:
        eps =  args['epochs']
    else:
        eps = 10
    if args['overwrite'] == 'True':
        overwrite = True
    else:
        overwrite = False


    if args['savedir']:
        savedir = args['savedir']
        if not os.path.exists(savedir):
            os.makedirs(savedir)
    else:
        savedir = './'
        # load json and create model
    if args['chunksize']:
        chunksize = args['chunksize']
    else:
        chunksize = 25

    if args['forest'] and args['hashes']:
        p = profiler.Profiler( hashes = args['hashes'], forest = args['forest'] , oma = True)
    else:
         p = profiler.Profiler( config_utils.datadir +'allnewlshforest.pkl', config_utils.datadir +'allhashes.h5',  oma = True)
    #shuffle
    df['HogFamA'] = df.HogA.map(hashutils.hogid2fam)
    df['HogFamB'] = df.HogB.map(hashutils.hogid2fam)

    df = df.sample(frac =1)
    msk = np.random.rand(len(df)) < 0.90

    #split
    traindf = df[msk]
    testdf = df[~msk]

    print( traindf)
    print( testdf)


    with open( config_utils.datadir + 'taxaIndex.pkl', 'rb')as taxain:
        taxaIndex = pickle.loads( taxain.read() )
    # 3 events, diff and union of matrows

    hogmat_size = 3 *2*  len(taxaIndex)

    if overwrite ==False & os.path.isfile(savedir + 'model.json') and  os.path.isfile(savedir + 'model.h5'):
        json_file = open(savedir + 'model.json', 'r')
        loaded_model_json = json_file.read()
        json_file.close()
        model = model_from_json(loaded_model_json)
        # load weights into new model
        model.load_weights(savedir + "model.h5")
        print("Loaded model from disk")
    else:
        layers = []
        layers.append(Dense( 10, input_dim= hogmat_size , activation='sigmoid' , use_bias=True))
        layers.append(Dense( 1, activation='softmax' , use_bias=True ))
        model = Sequential(layers)
    model.compile(loss='binary_crossentropy', optimizer='adam', metrics=['accuracy'])

    """

    tstart = t.time()
    testtrue = testdf.truth
        if overwrite == False & os.path.isfile(savedir + 'testset.pkl') :
            with open( savedir +  'testset.pkl' , 'rb') as testin:
                testdf = pickle.loads( testin.read())
            print(testdf)
        else:
            testset = list(set(list(testdf.HogFamA.unique()) + list(testdf.HogFamB.unique() ) ) )
            print(len(testset))
            #calculate the test set
            retdf = None
            for i in range( 0,len(testset),chunksize):
                print(i)
                retdict = p.retmat_mp(testset[i:i+chunksize])
                ret= pd.DataFrame.from_dict(retdict , orient= 'index')
                print(testset[i:i+chunksize])
                if i ==0 :
                    retdf = ret
                else:
                    retdf = pd.concat([retdf , ret])

            testdf = testdf.merge( retdf , left_on = 'HogFamA' , right_index = True , how= 'left')
            testdf = testdf.merge( retdf , left_on = 'HogFamB' , right_index = True , how= 'left')
            testdf = testdf.dropna(subset=['mat_y', 'mat_x'] , how = 'any')
            testdf['xtrain'] = testdf.apply( calculate_x , axis = 1)
            with open( savedir +  'testset.pkl' , 'wb') as testin:
                testin.write( pickle.dumps(testdf))

        X_test = np.vstack(testdf.xtrain)
        y_test = testdf.truth
    """
    tstart= t.time()
    gendata= p.retmat_mp(traindf, nworkers = 25, chunksize=50)
    #calculate testdf
    for X,y in gendata:
        print(X.shape)
        print(y.shape)
        metrics = model.train_on_batch(X, y )
        print(metrics)
        if t.time() - tstart > 200:
            # serialize model to
            model_json = model.to_json()
            with open(savedir + "model.json", "w") as json_file:
                json_file.write(model_json)
            # serialize weights to HDF5
            model.save_weights(savedir + "model.h5")
            print("Saved model to disk")
            tstart = t.time()
