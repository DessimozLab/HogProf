import numpy as np
from keras.models import *
from keras.optimizers import *
from keras.layers import *
from keras.metrics import *
from keras.regularizers import *
from keras.callbacks import *
import tensorflow as tf
import pickle
import collections
import itertools
from scipy.stats import bernoulli
import math
import random

import pandas as pd
import numpy as np
from qtlsearch.OBOParser import OBO
from goatools import obo_parser

import datetime
import gzip


omagaf =  '/home/cactuskid13/mntpt/OMA/latest/oma-go.txt'
obo = './go.obo'
unigaf = './goa_no_untrusted_iea.gaf.gz'

obofile = './go.obo'
obo = obo_parser.GODag(obofile)
qtl_go = OBO(obofile, store_as_int=True)
#redo only on OMA data w full ancestor nterms

def yeildBags(gaf, Yancestors= False):
    with open(gaf , 'r') as gafin:
        omaID = None
        lastID = None
        for l in gafin:
            if l[0] != '#':
                words = l.split()
                omaID,GO,evi,ref = words[0:4]
                if lastID is None:
                    lastID = omaID
                    if Yancestors == False:
                        GOset = set([int(GO.split(':')[1])])
                    else:
                        GOset = set(ancestors(GO) )

                    dat = { 'ID': omaID , 'GO':  GOset , 'REF' : [ref] , 'EVI' : [evi]}
                elif omaID != lastID:
                    yield dat
                    if Yancestors == False:
                        GOset = set([int(GO.split(':')[1])])
                    else:
                        GOset = ancestors(GO)

                    dat = { 'ID': omaID , 'GO': GOset , 'REF' : [ref] , 'EVI' : [evi]}
                else:
                    #todo: yeild ancestors of terms
                    if Yancestors == True:
                        dat['GO'] = dat['GO'].union( ancestors(GO) )
                    else:
                        dat['GO'] = dat['GO'].union(set([int(GO.split(':')[1])]))

                    dat['REF']+=[ref]
                    dat['EVI']+=[evi]
                lastID = omaID

def id2GO(intID):
    return 'GO:{:07d}'.format(intID)

def go2id(strID):

    if type( strID) is int:
        return strID
    try:
        return int(strID.split(':')[1])
    except:
        return int(strID)


def get_ancestors( golist , obo,  levels = 1 , verbose = True):
    #grab 4 levels up
    ancestors = []

    for l in range( levels):
        try:
            if l == 0 :
                a = [ obo[p.id].id for g in golist for p in obo[g].parents if g in obo and p.id in obo ]
            else:
                a = [ obo[p.id].id for g in a for p in obo[g].parents if g in obo and p.id in obo ]
            ancestors += a
        except:
            pass
    ancestors += golist
    return [ go2id(goterm) for goterm in set(ancestors) ]


def yeildBags_omagaf(gaf, obo, Yancestors= False):
    for omaid in gaf:
        yield { 'ID': omaid , 'GO':  get_ancestors(gaf[omaid] , obo ) }

def ancestors(term,verbose = True):
    try:
        numterm = int(term.split(':')[1])
    except :
        numterm = term
    try:
        ret = set(go.parents(numterm)).union(set([numterm]))
        return ret
    except:
        return set([numterm])

def buildGAF(gaf_file , universe= None):
    gaf_filtered = {}
    with open(gaf_file, mode='r') as gafin:

        for line in gafin:
            if line[0]!='#':
                words = line.split()
                if words[0] not in gaf_filtered:
                    gaf_filtered[words[0]]=set([words[1]])
                else:
                    gaf_filtered[words[0]].add(words[1])
    return gaf_filtered

def yeild_annotations_v2(gaf , obo):
    for prot in gaf:
        for a in set( get_ancestors( gaf[prot] , obo ) ):
            yield  a

def yeild_annotations(gaf, obo, verbose = True , Yancestors= False ):
     for entry in gaf:
        if Yancestors:
             for a in  get_ancestors( gaf[entry] , obo):
                 yield a
        else:
            for a in gaf[entry]:
                yield a

def makeGOdict(gaf , obo, sampling_factor= 1e-05 , Yancestors= True , uni = False):
    #return some descriptors of the dataset to be
    #used for embedding NN
    if uni == False:
        c = collections.Counter(yeild_annotations_v2(gaf, obo ))
    else:
        c = collections.Counter(yeild_annotations_uni(gaf))


    #count all the go terms in the OMA corpus
    nannot = sum(c.values())
    nterms = len(c.keys())
    #info = np.log(np.asarray(c.values())/nannot)
    #infocontent = dict(zip( c.keys(), list(info)))
    index = dict(zip(c.keys(), list(np.arange(nterms) )))
    reverse_index = dict( zip( index.values(), index.keys() ))
    freq = list(np.array(list(c.values()))/nannot)
    freq = [  (min(1, math.sqrt(word_frequency / sampling_factor) / (word_frequency / sampling_factor)  )  ) for word_frequency in freq ]
    sampling = dict(zip(index.values() , freq))
    return nterms , c , index , reverse_index , sampling

def prunesamples(samples , sampling , verbose = False):
    #remove samples in probabilistic way
    #select underepresented go terms in positive sampling
    #thin out pairs with overrepresented go terms in them
    if verbose == True:
        print([ sampling[ s[0] ] for s in samples ] )
        print( [ sampling[ s[1] ] for s in samples ] )
    ar1 = np.array([ random.uniform(0, 1) <= p  for p in [ sampling[ s[0] ] for s in samples ] ] , dtype = np.bool_ )
    ar2 = np.array([ random.uniform(0, 1) <= p  for p in [ sampling[ s[1] ] for s in samples ] ] , dtype = np.bool_ )
    select = np.bitwise_and(ar1,ar2)
    if np.sum(select)>0:
        samples = np.array(samples)[select,:]
        return samples
    else:
        return samples

def prunesamples_neg(samples , sampling ):
    #remove samples in probabilistic way
    #select overrepresented GO terms more often in negative sampling
    #ar1 = np.array([ random.uniform(0, 1) > p  for p in [ sampling[ s[0] ] for s in samples ] ] , dtype = np.bool_ )
    select = np.bitwise_and(ar1,ar2)
    samples = np.array(samples)[select,:]
    return samples

def yeildBags_uni(gaf, yeildancestors = True):
    with gzip.open(gaf,'rt') as gafin:
        ID = None
        lastID = None
        for l in gafin:
            if l[0] != '!':
                words = l.rstrip("\n").split("\t")

                ID,GO = (words[1], words[4])
                if lastID is None:
                    lastID = ID
                    dat = { 'ID': ID , 'GO': ancestors(GO) }
                elif ID != lastID:
                    yield dat
                    dat = { 'ID': ID , 'GO': ancestors(GO) }
                elif dat['GO']:
                    try:
                        dat['GO'] = dat['GO'].union( ancestors(GO) )
                    except TypeError:
                        pass
                else:
                    dat = { 'ID': ID , 'GO': ancestors(GO) }
                lastID = ID

def yeild_annotations_uni(gaf, verbose = False):
    with gzip.open(gaf,'rt') as gafin:
        for l in gafin:
            if l[0] != '!':
                #todo: yeild ancestors of terms
                #for term in retgoterms(l.split()[1]):
                #yield term
                try:
                    term = l.rstrip("\n").split("\t")[4]
                    for t in ancestors(term):
                        yield t
                except:
                    print(l)

def makesamples( gaf , obo , sampling, count ,  index , Yancestors = False , uni = False):
    #generator function to loop through gaf generating samples...
    if uni ==False:
        gafreader = yeildBags_omagaf(gaf, obo, Yancestors)
    else:
        gafreader = yeildBags_uni(gaf,   Yancestors)
    terms = list(index.keys())
    negatives = []
    print(sum(c.values()))
    pow = 1.5

    while len(negatives)< 10000000:
        neg1 = [ random.choice(terms) for i in range(100000) ]
        neg1 = [ n for n in neg1 if count[n]>0 and random.uniform(0, 1) < sampling[index[n]]**pow and count[n] > 50 ]
        negatives +=neg1
    #at least 100 annot in corpus
    thresh = 1

    infinite_gaf = itertools.cycle(gafreader)

    for i,dat in enumerate(infinite_gaf):
        try:
            #if i == 0:
            #    last =  [  index[s] for s in dat['GO'] if sampling[index[s]]**pow  < random.uniform(0,1) and count[index[s]] > thresh ]
            if len(dat['GO'])>1 and i > 0:
                samples = []
                maxiter = 100
                i = 0
                while len(samples) <1 and i < maxiter:
                    #favor less common words
                    samples =  [  index[s] for s in dat['GO'] if s in index and sampling[index[s]] > random.uniform(0,1)  and  count[s] > 20  ]
                    i += 1
                if i == maxiter:
                    samples =  [  index[s] for s in dat['GO'] if s in index ]

                posi = np.array([  [  c[0] ,c[1]  ]+ [1]  for c in itertools.combinations( samples , 2 )  if  c[0] != c[1] ] )
                nega = np.array([ [  random.choice(negatives) , random.choice(negatives) ] + [0]  for i in range(posi.shape[0])  ]  )

                samples =  np.vstack([posi,nega])
                if samples.shape[1]>1:
                    x1 = samples[:,0]
                    x2 = samples[:,1]
                    y = samples[:,2]
                    yield [x1,x2],y
            else:
                pass
        except ValueError:
            pass

Yancestors = True
readGAF = False

retrain = False
stamp = datetime.datetime.now().strftime('%m_%d_%Y_%H_%M_%S_%f')
min_annot= 50
gaf = buildGAF(omagaf)


if readGAF == True:
    #read the gaf to get term frequences
    nterms , c , index , reverse_index , sampling = makeGOdict( gaf , obo, sampling_factor= 1e-05 , Yancestors= Yancestors)
    print('done reading GAF')
    print('n terms:' + str(nterms))

    if Yancestors == True:
        with open( '../models/gafobects.pkl' , 'wb' ) as gafstats:
            gafstats.write(pickle.dumps([nterms , c , index , reverse_index , sampling]) )
    if Yancestors == False:
        with open( '../models/gafobects_noancestors.pkl' , 'wb' ) as gafstats:
            gafstats.write(pickle.dumps([nterms , c , index , reverse_index , sampling]) )

if Yancestors == True:
    with open( '../models/gafobects.pkl' , 'rb' ) as gafstats:
        nterms , c , index , reverse_index , sampling = pickle.loads(gafstats.read())
    modelfile = '../models/GO2vec'+stamp+'.h5'

if Yancestors == False:
    with open('../models/gafobects_noancestors.pkl' , 'rb' ) as gafstats:
        nterms , c , index , reverse_index , sampling = pickle.loads(gafstats.read())
    modelfile = '../models/GO2vec_noancestors'+stamp+'.h5'
print('train on n annot = ')

config = tf.ConfigProto()
config.gpu_options.per_process_gpu_memory_fraction= 0.95
K.set_session(tf.Session(config=config))

#filter the index
index = { t:index[t] for t in index if c[t] > min_annot }
nterms = len(index)

if retrain == False:
    #dimensionality of GO space
    vector_dim = 5
    #word2vec model to be trained
    input_target = Input((1,) , name='target_in')
    input_context = Input((1,) , name='context_in')

    embedding = Embedding(nterms, vector_dim, input_length=1, name='embedding' , embeddings_initializer='uniform', embeddings_regularizer= None,
     activity_regularizer=None, embeddings_constraint=None, mask_zero=False )

    target = embedding(input_target)
    target = Reshape((vector_dim, 1), name='target')(target)
    context = embedding(input_context)
    context = Reshape((vector_dim, 1) , name='context' )(context)
    similarity = dot([target, context], axes=0 , normalize = True )
    # now perform the dot product operation to get a similarity measure
    dot_product = dot([target, context] , axes=1)
    dot_product = Reshape((1,))(dot_product)

    # add the sigmoid output layer
    output = Dense(1, activation='sigmoid' , name = 'out')(dot_product)

    # create the primary training model
    #o = Adagrad(lr=0.001)
    #o = Adadelta(lr=1.0, rho=0.95)

    o = RMSprop(lr=0.025, rho=0.9)


    #o = Adagrad(lr=0.000075)

    model = Model(inputs=[input_target,input_context], outputs=[output])
    model.compile(loss='binary_crossentropy', optimizer=o , metrics = [ 'binary_accuracy'])
    embedder = Model( inputs=[input_target], output=target )
    validation_model = Model(input=[input_target, input_context], output=similarity)

    class SimilarityCallback:
        def run_sim(self):
            for i in range(valid_size):
                valid_word = reverse_dictionary[valid_examples[i]]
                top_k = 8  # number of nearest neighbors
                sim = self._get_sim(valid_examples[i])
                nearest = (-sim).argsort()[1:top_k + 1]
                log_str = 'Nearest to %s:' % valid_word
                for k in range(top_k):
                    close_word = reverse_dictionary[nearest[k]]
                    log_str = '%s %s,' % (log_str, close_word)
                print(log_str)

        @staticmethod
        def _get_sim(valid_word_idx):
            sim = np.zeros((vocab_size,))
            in_arr1 = np.zeros((1,))
            in_arr2 = np.zeros((1,))
            for i in range(vocab_size):
                in_arr1[0,] = valid_word_idx
                in_arr2[0,] = i
                out = validation_model.predict_on_batch([in_arr1, in_arr2])
                sim[i] = out
            return sim
    sim_cb = SimilarityCallback()
    ###modify this
    batchiter = 10000
    epochs = 100


if retrain == True:
    model = print('Load the model..')
    model = load_model(modelfile)
    #o = RMSprop(lr=0.0001, rho=0.9)
    o = Adadelta(lr=0.0001)
    #o = Nadam(lr=0.002, beta_1=0.9, beta_2=0.999)
    model.compile(loss='binary_crossentropy', optimizer=o , metrics = [ 'binary_accuracy'])

mc = ModelCheckpoint(modelfile, monitor = 'loss', mode = 'min', verbose = 1, save_best_only = False)
lr = ReduceLROnPlateau(monitor='loss', factor=0.5, patience= 20 , min_lr=0.000001 , verbose = 1)
tb = TensorBoard(log_dir='./logs', histogram_freq=0, batch_size=32, write_graph=True, write_grads=False, write_images=False, embeddings_freq=0, embeddings_layer_names=None, embeddings_metadata=None, embeddings_data=None, update_freq='epoch')

history = model.fit_generator(makesamples( gaf, obo ,sampling , c, index , Yancestors ), steps_per_epoch=10000, epochs=10000, verbose=1,
callbacks=[ mc , lr , tb  ], max_queue_size=10, workers=1, use_multiprocessing=False, shuffle=True, initial_epoch=0)

model.save(modelfile)
embedder.save(modelfile+'embedder')
