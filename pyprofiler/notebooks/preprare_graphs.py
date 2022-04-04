from pyprofiler.utils import hashutils
import rando
from pyprofiler.utils import config_utils
import pyprofiler.utils.goatools_utils as goa
import pyprofiler.utils.hashutils as hashutils
import numpy as np
import pyprofiler.profiler as profiler
import pandas as pd
import itertools
import redis
##to get the mapping of oma hogs to cogs to interactions in specific species I used dask distributed and a redis server
#you may need to get these up and running for you own cluster configuration before this notebook will work for you
import dask
import functools

import dendropy



import smallpars
import copy

import scipy
from dask import dataframe as dd
import pickle
from bloom_filter2 import BloomFilter
from sklearn.model_selection import train_test_split

#use the COG links to prepare phylo graphs using OMA and string data
#this is the input data for the phylo graphnet


#all string links 
string_links = '/scratch/dmoi/datasets/STRING/protein.links.full.v11.5.txt'

#coglinks with OMA fam IDS mapped. This is created in the notebook
prepared_coglinks = 'STRINGCOGS2OMAHOGS.csv'
stringPairs = pd.read_csv(prepared_coglinks)


#derive explicit profiles for our hogs of interest in string
calc_hogs_string = False
#calculate bloom filter for interaction presence between cogs in species
calc_filter = False
#use distributed computation to compile bloom filter
distributed = True



#add to existing data or overwrite
reload_dataset = False

#generate training data
traindata_gen = True
trainsample= 10000


#generate testing data
#create testing set using the testing set dataframe
testdata_gen = True
testsample= 2000



#open a profiler object
p = profiler.Profiler(lshforestpath = '/scratch/dmoi/datasets/all/newlshforest.pkl' , hashes_h5='/scratch/dmoi/datasets/birds/all/hashes.h5' , mat_path= None, oma = '/scratch/dmoi/datasets/OMA/apr2021/OmaServer.h5', tar= None , nsamples = 256 , mastertree = '/scratch/dmoi/datasets/birds/all_test_master_tree.pkl')

stringprofiles = {}
if calc_hogs_string == True:

	#lets load a compiled db containing the OMA root HOGs into a profiler oject 
	
	def grabHog(ID, verbose = True):
	    try:
	        entry = p.db_obj.entry_by_entry_nr(p.db_obj.id_resolver.resolve(ID))
	        return entry[4].decode() , entry
	    except:
	        return np.nan,np.nan
	


    print('profiles to calclulate',len(stringHOGs))
    for i,fam in enumerate(stringHOGs):
        if i % 100 ==0:
            print(i)
        try:
            prof = p.return_profile_OTF(fam)
            stringprofiles.update(prof)
        except:
            print('err',fam)
    with open('/scratch/dmoi/datasets/STRING/' + 'gold_standard_profiles.pkl' , 'wb') as profiles_out:
        profiles_out.write(pickle.dumps(stringprofiles))

else:
	with open('/scratch/dmoi/datasets/STRING/' + 'gold_standard_profiles.pkl' , 'rb' )as profiles_out:
    stringprofiles = pickle.loads(profiles_out.read())

string_df = pd.DataFrame.from_dict(stringprofiles , orient='index')

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


    if distributed == True:
        NCORE = 4
        print('deploying cluster')
        cluster = SLURMCluster(
            walltime='4:00:00',
            n_workers = NCORE,
            cores=NCORE,
            processes = NCORE,
            interface='ib0',
            memory="120GB",
            env_extra=[
            'source /scratch/dmoi/miniconda/etc/profile.d/conda.sh',
            'conda activate ML2'
            ],
            #scheduler_options={'interface': 'ens2f0' },
            #if gpu node
            scheduler_options={'interface': 'ens3f0' },
            #extra=["--lifetime", "3h55m", "--lifetime-stagger", "4m"]
        )
        print(cluster.job_script())
        print(cluster)
        cluster.scale(jobs = 100)
        print(cluster.dashboard_link)
        client = Client(cluster , timeout='450s' , set_as_default=True )
    else:
        cluster = LocalCluster()
        client = Client(cluster)

    link_df = dd.read_csv(string_links,  blocksize=75e6 , header = 0, sep = ' ')
    print(link_df)
    
	#compute bloom filters for protein pairs
	@dask.delayed
	def mapcogs(df ):
	    rdb = redis.Redis(host='10.202.12.174', port=6379, db=0)
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


    #b=BloomFilter(max_elements=10**8, error_rate=0.001 ,start_fresh = True)
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


with open('bloomfinal_big.pkl' , 'rb' ) as finalout:
    resfinal = pickle.loads(finalout.read()) 
print(resfinal)

#define a function that checks all the filters
def check_filters(element,filters):
    for f in filters:
        if element in f[0][0]:
            return True
    return False
bfilter = functools.partial(check_filters , filters= resfinal)

#sanity check
print(bfilter('COG1756_COG0088_4113'))
print(bfilter('crap'))

#open up the tax tree and convert to dendropy for use w fitch algo

dendrotree = dendropy.Tree.get(
        data=p.tree.write(format=3),
        rooting="force-rooted",
        suppress_internal_node_taxa=False,
        schema='newick')


#we're checking for interaction in a subset of species and propagating up
allowed_symbols =set([0,1,None])
transition_dict = { (c[0],c[1]):i for i,c in enumerate(itertools.permutations(allowed_symbols,2) ) }

#
def calc_interaction_on_taxonomy(cog1,cog2,treein , verbose = False):
    #set interaction states
    #look for interactions in bloom
    tree = copy.deepcopy(treein)
    for i,l in enumerate(tree.leaf_nodes()):
        l.event = {}
        l.scores = {}
        l.symbols = {}
        l.scores = { c:10**10 for c in allowed_symbols }
        if bfilter(cog1+'_'+cog2+'_'+l.taxon.label) or bfilter(cog2+'_'+cog1+'_'+l.taxon.label):
            if verbose == True:
                print(l.taxon.label)
            l.symbols = {1}
            l.scores[1] = 0
        else:
            l.symbols = {0}
            l.scores[0] = 0
        l.char = min(l.scores, key=l.scores.get)
    t = smallpars.calculate_small_parsimony(tree ,allowed_symbols, transition_dict)
    labels = np.array( [ n.char for n in t.nodes() ] )
    return  labels

#mapfam to matrow

string_fam_map = { f:i for i,f in enumerate(string_df.index)}
string_mat = np.vstack(string_df.mat)
print(string_mat.shape)
#train test split
Datasets = {}
for label,df,mapping,profilemat in [  ('string', stringPairs, string_fam_map ,string_mat )]: 
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

    
def get_sectors(tree, breakpt = 10):
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

def create_data_updown( tree, coglinkdf, profiles , taxindex , posi_percent = .5 ,  q = None , iolock= None,  verbose = False, loop= True ):
        #upward and downward connected phylo nodes
        connectmat_up, connectmat_down, connectmat_diag, template_features = tree2Single_sparse_graph_updown(tree)
        
        #sector based aggregation
        sectormat = get_sectors(tree, breakpt = 15 )
        num_total_nodes = template_features.shape[0]
        if iolock:
            with iolock:
                print('dataset size' , len(coglinkdf))
        else:
            print('dataset size' , len(coglinkdf))
        
        N = len(tree.nodes())
        Nsectors = sectormat.shape[1]
        allfams = list(set(coglinkdf.fam1.unique()).union( set(coglinkdf.fam2.unique() ) ))
        leafset = set([n.taxon.label for n in tree.leaf_nodes()])
        
        while True:
            toss = scipy.stats.bernoulli.rvs(posi_percent, loc=0, size=1, random_state=None)
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
                losses = [ n.name  for n in tp.traverse() if n.lost ]
                dupl = [ n.name  for n in tp.traverse() if n.dupl ]
                presence = [ n.name  for n in tp.traverse() if n.nbr_genes > 0   ]
                presences.append(presence)
                losses = [taxindex[n] for n in losses] 
                dupl = [taxindex[n] for n in dupl]
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
                    
                    #filter sectormat
                    subsectors = sectormat[matrows,:]
                    #use only connected sectors
                    i=np.where(subsectors.sum(axis=0)!=0)[1]
                    subsectors = subsectors[: , i]
                    Nsectors = subsectors.shape[1]
                    subsectors = sparse2pairs(subsectors)
                    #connect god node to sectors
                    
                    suboverview = scipy.sparse.lil_matrix((Nsectors , 2 ) )
                    suboverview[:,0] = 1
                    suboverview = sparse2pairs(suboverview)
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
                    data['sectornode'].x =torch.tensor(  np.zeros(( Nsectors ,sub_node_features.shape[1]))  ,  dtype=torch.double )
                    data['godnode'].x =torch.tensor(  np.zeros((1,sub_node_features.shape[1]))  ,  dtype=torch.double )

                    #up down fitch net
                    data['phylonodes_up', 'phylolink_up', 'phylonodes_up'].edge_index = torch.tensor(subconnect_up ,  dtype=torch.long )
                    data['phylonodes_down', 'phylolink_down', 'phylonodes_down'].edge_index = torch.tensor(subconnect_down ,  dtype=torch.long )             
                    data['phylonodes_up', 'phylolink_up_down', 'phylonodes_down'].edge_index = torch.tensor( subdiag ,  dtype=torch.long )
                    data['phylonodes_down', 'phylolink_down_up', 'phylonodes_up'].edge_index = torch.tensor( subdiag ,  dtype=torch.long )

                    #pooling connections
                    data['phylonodes_down', 'informs', 'sectornode'].edge_index = torch.tensor(subsectors,  dtype=torch.long )
                    data['phylonodes_up', 'informs', 'sectornode'].edge_index = torch.tensor(subsectors ,  dtype=torch.long )
                    
                    #final decision node
                    data['sectornode', 'informs', 'godnode'].edge_index = torch.tensor(suboverview,  dtype=torch.long )
                    
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

gen = create_data_updown( dendrotree, Datasets['string']['Train']  , stringprofiles , profile_mapper, .5  , verbose = False  )
print(next(gen))


#create training set using the generator on training samples
if traindata_gen == True:
    if reload_dataset == True:
        with open('trainingset.pkl' , 'rb')as trainout:
            samples = pickle.loads(trainout.read())
    else:
        print('newsamples')
        samples = []
    if traindata_gen == True:
        gen = create_data_updown( dendrotree, Datasets['string']['Train']  , stringprofiles , profile_mapper, .5  , verbose = False  )
        for i , data in enumerate( gen ):
            if i < 5:
                print(data)
                print(data['godnode'].y)
            if i % 100 ==0:
                print(i, len(samples))
                with open('trainingset.pkl' , 'wb')as trainout:
                    trainout.write(pickle.dumps(samples))
            samples.append(data)

if testdata_gen == True:
    if reload_dataset == True:
        with open('testgset.pkl' , 'rb')as trainout:
            samples = pickle.loads(trainout.read())
    else:
        print('newsamples')
        samples = []
    if testdata_gen == True:
        gen = create_data_updown( dendrotree, Datasets['string']['Test']  , stringprofiles , profile_mapper, .5  )
        for i , data in enumerate( gen ):
            if i % 100 ==0:
                print(i, len(samples))
                with open('testgset.pkl' , 'wb')as trainout:
                    trainout.write(pickle.dumps(samples))
            samples.append(data)

