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

linkfile = '/scratch/dmoi/datasets/STRING/protein.links.detailed.v11.5.txt'


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
      scheduler_options={'interface': 'ens2f0' },
      #extra=["--lifetime", "3h55m", "--lifetime-stagger", "4m"]
  )
  print(cluster.job_script())

else:
  cluster = LocalCluster()
  client = Client(cluster)

if distributed == True:
  print(cluster)
  cluster.scale(jobs = 100)
  print(cluster.dashboard_link)
  client = Client(cluster , timeout='450s' , set_as_default=True )

#find which species each of the cogs has an interaction in

#link_df = dd.read_csv('/scratch/dmoi/datasets/STRING/protein.physical.links.detailed.v11.5.txt', blocksize=100e6 , header = 0, sep = ' ')
link_df = dd.read_csv(linkfile,  blocksize=100e6 , header = 0, sep = ' ')
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
   return   b , len(coglinks)

@dask.delayed
def sumfilter(f1,f2, total ):
   if type( f1 ) == tuple:
     f1 = f1[0]
   if type( f2 ) == tuple:
     f2 = f2[0]
   f3 = f1.__ior__(f2)

   return f3 , total

partitions  = link_df.to_delayed()
print('map cogs')
res1 = [ mapcogs(p) for p in partitions ]
print('done')
print('make filters')
res2 = [ return_filter(p) for p in res1 ] 
totalfilter = res2
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
print('done')
resfinal = dask.compute(totalfilter)
if calc_filter == True:
 with open('bloomfinal_big.pkl' , 'wb' ) as finalout:
     finalout.write(pickle.dumps(resfinal))