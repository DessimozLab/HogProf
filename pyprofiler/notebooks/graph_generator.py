#!/usr/bin/env python
# coding: utf-8

# In[1]:


import sys
sys.path.append('../..')
#sys.path.append( '/home/cactuskid13/miniconda3/pkgs/')
print(sys.path)


# In[2]:


import torch


# In[3]:


from pyprofiler.utils import hashutils
import ete3
import random
from pyprofiler.utils import config_utils
import pyprofiler.utils.goatools_utils as goa
import pyprofiler.utils.hashutils as hashutils
import seaborn as sns
from matplotlib import pyplot as plt
import numpy as np
import pyprofiler.profiler as profiler
import pandas as pd
import itertools
import redis
##to get the mapping of oma hogs to cogs to interactions in specific species I used dask distributed and a redis server
#you may need to get these up and running for you own cluster configuration before this notebook will work for you
import dask
import warnings


from sklearn.metrics import precision_recall_curve
import matplotlib.pyplot as plt
from sklearn.metrics import roc_curve
from sklearn.metrics import auc


import scipy
from dask import dataframe as dd
import pickle
from bloom_filter2 import BloomFilter
from sklearn.model_selection import train_test_split


# In[4]:


#lets load a compiled db containing the OMA root HOGs into a profiler oject 
p = profiler.Profiler(lshforestpath = '/scratch/dmoi/datasets/all/newlshforest.pkl' , hashes_h5='/scratch/dmoi/datasets/birds/all/hashes.h5' , mat_path= None, oma = '/scratch/dmoi/datasets/OMA/apr2021/OmaServer.h5', tar= None , nsamples = 256 , mastertree = '/scratch/dmoi/datasets/birds/all_test_master_tree.pkl')


# In[5]:


def grabHog(ID, verbose = True):
    try:
        entry = p.db_obj.entry_by_entry_nr(p.db_obj.id_resolver.resolve(ID))
        return entry[4].decode() , entry
    except:
        return np.nan,np.nan
#map to OMA HOGs


# In[6]:

filter_coglinks = False

if filter_coglinks == True:
    from collections import Counter
    coglink_df = dd.read_csv('/scratch/dmoi/datasets/STRING/COG.links.detailed.v11.5.txt', blocksize=25e6 , header = 0, sep = ' ')
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


# In[14]:
#map the interacting cogs to the proteins
compute_grabcogs = False
if compute_grabcogs == True:
    grabcogs = set( list(coglink_df.group1.unique()) + list(coglink_df.group2.unique()) )
    grabcogs= list(grabcogs)
    COGmapings_df = dd.read_csv('/scratch/dmoi/datasets/STRING/COG.mappings.v11.5.txt', blocksize=25e6 , header = 0, sep = '\t')
    COGmapings_df = COGmapings_df.set_index('orthologous_group')
    COGmapings_df.astype(str)
    COGmapings_df['##protein'].map( lambda x : x.strip() )
    COGmapings_df['species'] = COGmapings_df['##protein'].map( lambda x : x.split('.')[0] )
    COGmapings_df['COG'] = COGmapings_df.index
    COGmapings_df = COGmapings_df.loc[grabcogs]
    COGmapings_df = COGmapings_df.compute()
    print(COGmapings_df.head())


# In[15]:


#only take the proteins in our cogs of interest
if compute_grabcogs == True:
    grabprots =list(COGmapings_df['##protein'].unique())
    print(len(grabprots))
    with open('/scratch/dmoi/datasets/STRING/COG.links.detailed.v11.5.txt' + '.grabcogs.txt', 'w') as protsout:
        protsout.write(''.join([ p + '\n' for p in grabcogs ]) )
    with open('/scratch/dmoi/datasets/STRING/COG.mappings.v11.5.txt' + '.grabprots.txt' , 'w') as protsout:
        protsout.write(''.join([ p + '\n' for p in grabprots ]) )
else:
    with open('/scratch/dmoi/datasets/STRING/COG.links.detailed.v11.5.txt' + '.grabcogs.txt', 'r') as protsout:
        grabcogs = [ cog for cog in protsout.readlines()]
    with open('/scratch/dmoi/datasets/STRING/COG.mappings.v11.5.txt' + '.grabprots.txt' , 'r') as protsout:
        grabprots = [ prot for prot in protsout.readlines()]


# In[16]:


calc_mappers = False
if calc_mappers == True:
    rdb = redis.Redis(host='10.202.12.174', port=6379, db=0)

    count = 0
    for i,r in COGmapings_df.iterrows():
        rdb.set(r['##protein'], i)
        count+=1
        if count < 10:
            print(i+'\n',r)
        if count%1000000==0:
            print(count/len(COGmapings_df))


# In[17]:


maphogs = False
if maphogs == True:
    #mapping each string cog to an oma hog by selecting a member of the cog
    rdb = redis.Redis(host='10.202.12.174', port=6379, db=0)
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


# In[18]:


print(len(hogmap))
for i, key in enumerate(hogmap):
    if i < 10:
        print(key, hogmap[key])


# In[19]:


#add the HOGs to the COGdf
#grab the corresponding profiles
compile_final_cogdf = False
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


# In[20]:
stringPairs = coglink_df
# In[21]:
with open('/scratch/dmoi/datasets/STRING/' + 'gold_standard_profiles.pkl' , 'rb' )as profiles_out:
    stringprofiles = pickle.loads(profiles_out.read())

string_df = pd.DataFrame.from_dict(stringprofiles , orient='index')

#make the profiles for this small set of HOGs
for i, key in enumerate(stringprofiles):
    if i < 10:
        print(key,stringprofiles[key])


with open('bloomfinal_big.pkl' , 'rb' ) as finalout:
    resfinal = pickle.loads(finalout.read()) 
print(resfinal)


# In[35]:


def check_filters(element,filters):
    for f in filters:
        if element in f[0][0]:
            return True
    return False
import functools
bfilter = functools.partial(check_filters , filters= resfinal)
print(bfilter('COG1756_COG0088_4113'))
print(bfilter('crap'))


taxmapper = 'string2oma_specmap.pkl'
calc_taxmapper = True
if calc_taxmapper == True:
     with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        #compare lineages to get the closest
        from ete3 import NCBITaxa
        ncbi = NCBITaxa()
        #map profiler leaves to closest leaf in string
        strings_species = '/scratch/dmoi/datasets/STRING/species.v11.5.txt'
        string_specdata = pd.read_table(strings_species)
        print(string_specdata)
        stringset = set([ str(tax) for tax in list(string_specdata['#taxon_id']) ]) 
        omaset =  set([ spec.name  for spec in p.tree.get_leaves()])
        print('omaset',len(omaset))
        print('stringset',len(stringset))
        shared_leaves = omaset.intersection(stringset)
        missing_leaves = omaset - stringset
        #these are oma leaves not in string...we can map to the closest string species.
        string_missing = stringset - omaset
        omalineages = {tax: set(ncbi.get_lineage(int(tax))) for tax in missing_leaves}
        stringlineages = {tax: set(ncbi.get_lineage(int(tax))) for tax in string_missing}
        string2oma={}
        print('done lineages')
        for i,tax in enumerate(stringlineages):
            #find the closest
            shared = { tax_oma: len(stringlineages[tax].intersection(omalineages[tax_oma]))/len(omalineages[tax_oma]) for tax_oma in omalineages } 

            string2oma[tax] = max(shared, key=shared.get)
            if i %1000 == 0:
                print(i/ len(stringlineages) )
        with open(taxmapper, 'wb') as taxmapper:
            taxmapper.write(pickle.dumps(string2oma))
else:
    with open(taxmapper, 'rb') as taxmapper:
        string2oma= pickle.loads(taxmapper.read())



#test out to find the species for a cog pair
species = [ spec.name  for spec in p.tree.get_leaves()] + list(string2oma.keys())
species = set(species)


import dendropy
taxnwk = '/scratch/dmoi/datasets/birds/all_test_master_tree.nwk'
with open( 'taxtree.nwk' , 'w') as treeout:
    treeout.write(p.tree.write())
dendrotree = dendropy.Tree.get(
        data=p.tree.write(format=3),
        rooting="force-rooted",
        suppress_internal_node_taxa=False,
        schema='newick')


#setup the internal nodes for the fitch algo
for i,l in enumerate(dendrotree.nodes()):
    l.event = {}
    l.scores = {}
    l.symbols = None
    l.char= None
    l.matrow = i


import smallpars
import copy

#we're checking for interaction in a subset of species and propagating up
allowed_symbols =set([0,1,None])
transition_dict = { (c[0],c[1]):i for i,c in enumerate(itertools.permutations(allowed_symbols,2) ) }
def calc_interaction_on_taxonomy(cog1,cog2,treein ,species_set = species, string2oma= string2oma, verbose = False):
    #set interaction states
    #look for interactions in bloom
    coglink1 = cog1+'_'+cog2 +'_'
    coglinks_species = [ coglink1+spec for spec in species_set ]
    coglink2 = cog2+'_'+cog1+'_'
    coglinks_species += [ coglink2+spec for spec in species_set ]
    
    checklinks = [ bfilter( spec) for spec in coglinks_species ]
    
    species_set = [ s for s,c in list(zip(species_set,checklinks))  if c== True  ]
    if verbose == True:
        print(cog1,cog2)
        print('string entries:',len(species_set))
    species_set = set([ s if s not in string2oma else string2oma[s] for s in species_set  ])
    tree = copy.deepcopy(treein)
    for i,l in enumerate(tree.leaf_nodes()):
        l.event = {}
        l.scores = {}
        l.symbols = {}
        l.scores = { c:10**10 for c in allowed_symbols }
        if l.taxon.label in species_set:
            l.symbols = {1}
            l.scores[1] = 0
        else:
            l.symbols = {0}
            l.scores[0] = 0
        l.char = min(l.scores, key=l.scores.get)
    t = smallpars.calculate_small_parsimony(tree ,allowed_symbols, transition_dict)
    labels = np.array( [ n.char for n in t.nodes() ] )
    return  labels

string_fam_map = { f:i for i,f in enumerate(string_df.index)}
string_mat = np.vstack(string_df.mat)
print(string_mat.shape)

#train test split
Datasets = {}
for label,df,mapping,profilemat in [  ('string', stringPairs, string_fam_map ,string_mat ) ]: 
    keys = set(mapping.keys())
    entry1 = [ f in keys for f in df.fam1]
    df = df.iloc[entry1]
    entry2 = [ f in keys for f in df.fam2]
    df = df.iloc[entry2]
    msk = np.random.rand(len(df)) < 0.8
    df_train = df.iloc[msk]
    df_test = df.iloc[~msk]
    Datasets[label]={'Train':df_train,'Test':df_test , 'mapping': mapping , 'mat':profilemat }
    print(label)
    print('train',len(df_train))
    print('test',len(df_test))


# In[47]:


import torch.nn.functional as F
import torch_geometric.transforms as T
from torch_geometric.nn import ChebConv
from torch_geometric.nn import  to_hetero
from torch_geometric.data import Data
from torch_geometric.loader import DataLoader
from torch_geometric.data import HeteroData ,InMemoryDataset
import copy
import time
import pickle
import random


# In[54]:


nodes = set([ n.name for n in p.tree.traverse() ])
dendrotree_nodes = set([str(n.taxon.label) if n.taxon else '-1' for n in dendrotree.nodes()] )


# In[55]:


print(len(nodes))
print(len(dendrotree_nodes))
print(len(nodes.intersection(dendrotree_nodes)))


# In[56]:


profile_mapper = { n.name:i for i,n in enumerate(p.tree.traverse()) }


# In[57]:


profile_mapper = { (n.taxon.label if n.taxon else '-1'):n.matrow for n in  dendrotree.nodes() }

def tree2Single_sparse_graph_updown(tree):
    N = len(tree.nodes())
    #mimic the fitch algo
    #propagate up and down in separate graphs
    index_up = np.vstack([ [n.matrow, c.matrow ] for n in tree.nodes() for c in n.child_nodes()])
    index_down = np.vstack([ [c.matrow, n.matrow ] for n in tree.nodes() for c in n.child_nodes()])
    connectmat_up = scipy.sparse.lil_matrix(( N ,  N ) )
    connectmat_down = scipy.sparse.lil_matrix(( N ,  N ) )
    connectmat_up[index_up[:,0],index_up[:,1]] = 1 
    connectmat_down[index_down[:,0],index_down[:,1]] = 1 
    diag = [[n,n] for n in range(N)]
    connectmat_diag=scipy.sparse.lil_matrix(( N ,  N ) )
    connectmat_diag[diag,diag] = 1 
    ntime = np.array([ n.distance_from_root() for n in tree.nodes()])
    mtime = np.amax(ntime)
    ntime/=mtime
    levels = np.array([ n.level() for n in tree.nodes() ] , dtype='double')
    levels /= np.amax(levels)
    Norm_nchild= np.array( [ len(n.child_nodes()) for n in tree.nodes() ] ,dtype='double' )
    mchild =np.amax(Norm_nchild)
    Norm_nchild/=mchild 
    Norm_nsister= np.array( [ len(n.sister_nodes()) for n in tree.nodes() ] ,dtype='double' )
    msis =np.amax(Norm_nsister)
    Norm_nsister/=msis    
    template_features = np.stack([ntime ,  Norm_nchild , Norm_nsister ]).T    
    return connectmat_up, connectmat_down, connectmat_diag, template_features


# In[59]:



def getmrca(treein,taxset):
    
    tree = copy.deepcopy(treein)
    n = tree.mrca(taxon_labels=taxset)
    if n is not None:
        subtree = dendropy.Tree(seed_node=n)
        taxset = set([ t.taxon.label if t.taxon else '-1'  for t in subtree.nodes()])
        matrows = [ t.matrow for t in treein.nodes() if t.taxon and t.taxon.label in taxset]
        return taxset, matrows , n
    else: 
        return None, None, None



def sparse2pairs(sparsemat, matrows = None):
    if matrows :
        sparsemat = sparsemat[matrows,:]
        sparsemat = sparsemat[:,matrows]
    sparsemat = scipy.sparse.find(sparsemat)
    return np.vstack([sparsemat[0],sparsemat[1]])


def process_node_down(node, sector = 0, breakpt = 10 , total = 0 ):
    node.sector = sector
    if sector == 0 :
        global count
        count = 0
    total += len(node.child_nodes())
    for i,child in enumerate(node.child_nodes()):
        if total > breakpt:
            if len(child.child_nodes())>0:
                #new sector w new total
                count+=1
                process_node_down(child, count , total = 0 , breakpt = breakpt)
            else:
                #leaf
                process_node_down(child, count , total = 0 , breakpt = breakpt)
        else:
            process_node_down(child, count , total = total , breakpt = breakpt)

    
def get_sectors(tree, breakpt = 20):
    process_node_down( tree.seed_node , sector = 0, breakpt = breakpt )
    row = [n.matrow for n in tree.nodes()]
    col = [n.sector for n in tree.nodes()]
    data = np.ones((len(row)))
    sectormat = scipy.sparse.csc_matrix( (data,(row,col)) )
    return sectormat

for i,l in enumerate(dendrotree.nodes()):
    l.sum_lengths = None
for i,l in enumerate(dendrotree.leaf_nodes()):
    l.sum_lengths = 1


def create_data_updown_nosectors( tree, coglinkdf, profiles , taxindex , posi_percent = .5 ,  q = None , iolock= None,  verbose = False, loop= True ):
        #upward and downward connected phylo nodes
        connectmat_up, connectmat_down, connectmat_diag, template_features = tree2Single_sparse_graph_updown(tree)
        N = len(tree.nodes())
        allfams = list(set(coglinkdf.fam1.unique()).union( set(coglinkdf.fam2.unique() ) ))
        leafset = set([n.taxon.label for n in tree.leaf_nodes()])
        while True:
            toss = scipy.stats.bernoulli.rvs(posi_percent, loc=0, size=1, random_state=42)
            if verbose == True:
                print('posi/nega',toss)
            if toss == 0:
                fam1 = random.choice(allfams)
                fam2 = fam1
                while fam1 == fam2:
                    fam2 = random.choice(allfams)
                labels = np.zeros((template_features.shape[0],))
            else:
                #positive sample
                dfline = coglinkdf.sample(n=1, random_state = random.randint(1,1000)).iloc[0]
                cog1= str(dfline.group1).replace("b",'').replace("'",'').strip()
                cog2= str(dfline.group2).replace("b",'').replace("'",'').strip()
                fam1 = dfline.fam1
                fam2 = dfline.fam2
                labels = calc_interaction_on_taxonomy(cog1,cog2,tree)
            nodefeatures = []
            presences= []
            #find profile nameset
            for i,tp in enumerate([profiles[fam1]['tree'], profiles[fam2]['tree']]):    
                profilefeatures = np.zeros((template_features.shape[0],3) )
                #find on which nodes the events happened
                losses = [ taxindex[n.name]  for n in tp.traverse() if n.lost ]
                dupl = [ taxindex[n.name]  for n in tp.traverse() if n.dupl ]
                presence = [ n.name  for n in tp.traverse() if n.nbr_genes > 0   ]
                presences.append(presence)
                presence = [taxindex[n] for n in presence]
                profilefeatures[losses, 0] = 1
                profilefeatures[dupl, 1] = 1
                profilefeatures[presence, 2] = 1
                nodefeatures.append(profilefeatures)
            nodeset = set(presences[0]).union(set(presences[1]))
            if len(nodeset)> 10:
                skip = False
                try:
                    taxset,matrows,n = getmrca(tree,leafset.intersection(nodeset))
                    if taxset == None:
                        skip = True
                except ValueError:
                    #no species overlap
                    skip = True

                if skip == False:
                    
                    #pare down labels
                    labels = labels[matrows]
                    
                    if verbose == True:
                        print('features',nodefeatures)
                        print( 'labels' , labels)
                    neglabels = np.ones(labels.shape)
                    neglabels = neglabels - labels
                    labels = np.vstack([labels,neglabels]).T

                    overview = scipy.sparse.lil_matrix( (len(matrows) , 2 ) )
                    overview[:,0] = 1

                    overview_rev = sparse2pairs(overview.T)
                    overview = sparse2pairs(overview)

                    #phylonode connections
                    subconnect_up = sparse2pairs(connectmat_up, matrows)
                    subconnect_down = sparse2pairs(connectmat_down, matrows)
                    subdiag = sparse2pairs(connectmat_diag, matrows)
                    
                    
                    #profile features
                    nodefeatures=np.hstack(nodefeatures)
                    
                    sub_template_features= template_features[matrows,:]
                    sub_node_features= nodefeatures[matrows,:]
                    sub_node_features = np.hstack([sub_template_features , sub_node_features])
                    godlabel = np.ones((1,1))*toss
                    godlabel = np.hstack([np.ones((1,1))-godlabel, godlabel])
                    
                    
                    data = HeteroData()   
                    
                    data['phylonodes_up'].x = torch.tensor(sub_node_features, dtype=torch.double )
                    data['phylonodes_down'].x = torch.tensor(sub_node_features, dtype=torch.double )
                    data['godnode'].x =torch.tensor(  np.zeros((1,1))  ,  dtype=torch.double )

                    #up down fitch net
                    data['phylonodes_up', 'phylolink_up', 'phylonodes_up'].edge_index = torch.tensor(subconnect_up ,  dtype=torch.long )
                    data['phylonodes_down', 'phylolink_down', 'phylonodes_down'].edge_index = torch.tensor(subconnect_down ,  dtype=torch.long )             
                    data['phylonodes_up', 'phylolink_up_down', 'phylonodes_down'].edge_index = torch.tensor( subdiag ,  dtype=torch.long )
                    data['phylonodes_down', 'phylolink_down_up', 'phylonodes_up'].edge_index = torch.tensor( subdiag ,  dtype=torch.long )

                    #pooling connections
                    data['phylonodes_down', 'informs', 'godnode'].edge_index = torch.tensor(overview ,  dtype=torch.long )
                    data['phylonodes_up', 'informs', 'godnode'].edge_index = torch.tensor(overview ,  dtype=torch.long )
                    
                    #pooling connections
                    data['godnode',  'informs','phylonodes_down' ].edge_index = torch.tensor(overview_rev,  dtype=torch.long )
                    data['godnode',  'informs','phylonodes_up'].edge_index = torch.tensor(overview_rev ,  dtype=torch.long )
                    
                    #add 2 label classes                
                    data['phylonodes_up'].y =torch.tensor(labels  ,  dtype=torch.long )
                    data['phylonodes_down'].y =torch.tensor(labels  ,  dtype=torch.long )
                    data['godnode'].y =torch.tensor( godlabel  ,  dtype=torch.long )

                    data = T.AddSelfLoops()(data)
                    data = T.NormalizeFeatures()(data)
                    if q:
                        q.put(data)
                    else:
                        yield data


gen = create_data_updown_nosectors( dendrotree, Datasets['string']['Train']  , stringprofiles , profile_mapper, .5  , verbose = False  )
print(next(gen))

# In[ ]:

traindata_gen = True
trainsample= 10000



#create training set using the generator on training samples
reload = True
if traindata_gen == True:
    if reload == True:
        with open('trainingset_nosectors.pkl' , 'rb')as trainout:
            samples = pickle.loads(trainout.read())
    else:
        print('newsamples')
        samples = []
    if traindata_gen == True:
        gen = create_data_updown_nosectors( dendrotree, Datasets['string']['Train']  , stringprofiles , profile_mapper, .5  , verbose = False  )
        for i , data in enumerate( gen ):
            if i < 5:
                print(data)
                print(data['godnode'].y)
            if i % 100 ==0:
                print(i, len(samples))
                with open('trainingset_nosectors.pkl' , 'wb')as trainout:
                    trainout.write(pickle.dumps(samples))
            samples.append(data)
            if i > trainsample:
                break

'''
#create testing set using the testing set dataframe
testdata_gen = True
testsample= 5000
NCORE = 20
reload = True
if testdata_gen == True:
    if reload == True:
        with open('testgset_nosectors.pkl' , 'rb')as trainout:
            samples = pickle.loads(trainout.read())
    else:
        print('newsamples')
        samples = []
    if testdata_gen == True:
        gen = create_data_updown_nosectors( dendrotree, Datasets['string']['Test']  , stringprofiles , profile_mapper, .5  )
        for i , data in enumerate( gen ):
            if i % 100 ==0:
                print(i, len(samples))
                with open('testgset_nosectors.pkl' , 'wb')as trainout:
                    trainout.write(pickle.dumps(samples))
            samples.append(data)

            if i > testsample:
                break

'''