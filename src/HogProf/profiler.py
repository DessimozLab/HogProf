from pyoma.browser import db
import pickle
import pandas as pd
import h5py
import random
from tables import *
import numpy as np
import random
import ete3
#from validation import validation_semantic_similarity
from .utils import hashutils , pyhamutils , files_utils
from time import time
import multiprocessing as mp
import functools
import numpy as np
import time
import gc
import logging
from pyoma.browser import db

np.random.seed(0)
random.seed(0)

class Profiler:

	"""
	A profiler object allows the user to query the LSH with HOGs and get a list of result HOGs back

	"""
	def __init__(self,lshforestpath = None, hashes_h5=None, mat_path= None, oma = False , nsamples = 256 , mastertree = None ):
		"""
		The Profiler class initializes a profiler object for querying the LSH with HOGs and returning a list of result HOGs.

		Attributes:
		lshobj (object): LSH object for querying.
		hashes_h5 (h5py.File): H5 file containing HOGs.
		nsamples (int): Number of samples to use.
		tree (ete3.Tree): Master tree used for generating taxa index.
		tree_string (str): String representation of the master tree.
		taxaIndex (dict): Dictionary mapping taxa names to their indices in the master tree.
		ReverseTaxaIndex (dict): Dictionary mapping indices in the master tree to their corresponding taxa names.
		db_obj (db.Database): OMA database object.
		treeweights (dict): Dictionary containing the tree weight for each taxon.
		READ_ORTHO (callable): Function for reading orthoxml files from OMA.
		HAM_PIPELINE (callable): Function for generating the Annotated tree from a row.
		HASH_PIPELINE (callable): Function for generating the hash from a row.

		Parameters:
		lshforestpath (str, optional): Path to the pickled LSH forest object.
		hashes_h5 (str, optional): Path to the H5 file containing HOGs.
		mat_path (str, optional): Path to the matrix file containing HOGs.
		oma (str, optional): Path to the OMA database.
		tar (str, optional): Path to the tar archive.
		nsamples (int, optional): Number of samples to use. Defaults to 256.
		mastertree (str, optional): Path to the master tree file.
		"""

		print('loading lsh')
		with open(lshforestpath, 'rb') as lshpickle:
			self.lshobj = pickle.loads(lshpickle.read())
			print('indexing lsh')
			self.lshobj.index()

		self.hashes_h5 = h5py.File(hashes_h5, mode='r')
		self.nsamples = nsamples

		if mastertree.split('.')[-1] == 'pkl':
				with open( mastertree , 'rb') as pklin:
					self.tree = pickle.loads(pklin.read())
					self.tree_string = self.tree.write(format=1)
		elif mastertree.split('.')[-1] == 'nwk':
			self.tree = ete3.Tree(mastertree,format=1)
			self.tree_string = self.tree.write(format=1)
		
		else:
			raise Exception( 'please provide a pickled ete3 tree or a newick file' )
		self.taxaIndex, self.ReverseTaxaIndex = files_utils.generate_taxa_index(self.tree)
			
		if oma:
			h5_oma = open_file(oma, mode="r")
			self.db_obj = db.Database(h5_oma)
			self.treeweights = hashutils.generate_treeweights(self.tree , self.taxaIndex , None, None )
			self.READ_ORTHO = functools.partial(pyhamutils.get_orthoxml_oma	, db_obj=self.db_obj)
			self.HAM_PIPELINE = functools.partial(pyhamutils.get_ham_treemap_from_row, tree=self.tree_string )
			self.HASH_PIPELINE = functools.partial(hashutils.row2hash , taxaIndex=self.taxaIndex  , treeweights=self.treeweights , wmg=None )

		print('DONE')

	def hogid2fam(self, hog_entry):
		if type(hog_entry )== int:
			return hog_entry
		else:
			try:
				return int(self.db_obj.hog_family( entry=hog_entry ) )
			except:
				return np.nan

	def return_profile_OTF(self, fam):
		"""
		Returns profiles as binary vectors for use with optimisation pipelines
		"""
		if type(fam) is str:
			fam = self.hogid2fam(fam)
		ortho_fam = self.READ_ORTHO(fam)
		if ortho_fam:
			tp = self.HAM_PIPELINE([fam, ortho_fam])

			losses = [ self.taxaIndex[n.name]  for n in tp.traverse() if n.lost and n.name in self.taxaIndex  ]
			dupl = [ self.taxaIndex[n.name]  for n in tp.traverse() if n.dupl  and n.name in self.taxaIndex  ]
			presence = [ self.taxaIndex[n.name]  for n in tp.traverse() if n.nbr_genes > 0  and n.name in self.taxaIndex  ]

			indices = dict(zip (['presence', 'loss', 'dup'],[presence,losses,dupl] ) )
			hog_matrix_raw = np.zeros((1, 3*len(self.taxaIndex)))
			for i,event in enumerate(indices):
				if len(indices[event])>0:
					taxindex = np.asarray(indices[event])
					hogindex = np.asarray(indices[event])+i*len(self.taxaIndex)
					hog_matrix_raw[:,hogindex] = 1
			return {fam:{ 'mat':hog_matrix_raw, 'tree':tp} }
		else:
			return{ fam: { 'mat':None , 'tree':None }}


	def return_profile_complements(self, fam):
		"""
		Returns profiles for each loss to search for complementary hogs
		"""
		if type(fam) is str:
			fam = self.hogid2fam(fam)
		ortho_fam = self.READ_ORTHO(fam)
		tp = self.HAM_PIPELINE([fam, ortho_fam])

		losses = set([ n.name  for n in tp.traverse() if n.lost and n.name in self.taxaIndex  ])
		#these are the roots of the fams we are looking for
		#we just assume no duplications or losses from this point

		ancestral_nodes = ([ n for n in profiler.tree.traverse() if n.name in losses])
		losses=[]
		dupl=[]
		complements={ n.name+'_loss' : [] }

		indices = dict(zip (['presence', 'loss', 'dup'],[presence,losses,dupl] ) )

		hog_matrix_raw = np.zeros((1, 3*len(self.taxaIndex)))
		for i,event in enumerate(indices):
			if len(indices[event])>0:
				taxindex = np.asarray(indices[event])
				hogindex = np.asarray(indices[event])+i*len(self.taxaIndex)
				hog_matrix_raw[:,hogindex] = 1
		
		return {fam:{ 'mat':hog_matrix_raw, 'hash':tp} }

	def worker( self,i, inq, retq ):
		"""
		this worker function is for parallelization of generation of binary vector for use with optimisation pipelines

		"""
		print('worker start'+str(i))
		while True:
			input = inq.get()
			if input is None:
				break
			else:
				fam,ortho_fam = input
				tp = self.HAM_PIPELINE([fam, ortho_fam])
				losses = [ self.taxaIndex[n.name]  for n in tp.traverse() if n.lost and n.name in self.taxaIndex  ]
				dupl = [ self.taxaIndex[n.name]  for n in tp.traverse() if n.dupl  and n.name in self.taxaIndex  ]
				presence = [ self.taxaIndex[n.name]  for n in tp.traverse() if n.nbr_genes > 0  and n.name in self.taxaIndex  ]
				indices = dict(zip (['presence', 'loss', 'dup'],[presence,losses,dupl] ) )
				hog_matrix_raw = np.zeros((1, 3*len(self.taxaIndex)))
				for i,event in enumerate(indices):
					if len(indices[event])>0:
						taxindex = np.asarray(indices[event])
						hogindex = np.asarray(indices[event])+i*len(self.taxaIndex)
						hog_matrix_raw[:,hogindex] = 1
				retq.put({fam:{ 'mat':hog_matrix_raw, 'tree':tp} })


	def retmat_mp(self, traindf , nworkers = 25, chunksize=50  ):
		"""
		function used to create training matrix with pairs of hogs. calculate_x will return the intersetcion of
		two binary vectors generated by pyham
		"""
		#fams = [ hashutils.hogid2fam(fam) for fam in fams ]
		def calculate_x(row):
			mat_x1 = row.mat_x
			mat_x2 = row.mat_y
			ret1 = np.zeros(mat_x1.shape)
			ret2 = np.zeros(mat_x2.shape)
			#diff = mat_x1 - mat_x2
			matsum = mat_x1 + mat_x2
			#ret1[np.where(diff != 0 ) ] = -1
			ret2[ np.where(matsum == 2 ) ] = 1
			return list(ret2)
		retq= mp.Queue(-1)
		inq= mp.Queue(-1)
		processes = {}
		mp.log_to_stderr()
		logger = mp.get_logger()
		logger.setLevel(logging.INFO)

		for i in range(nworkers):
			processes[i] = {'time':time.time() , 'process': mp.Process( target = self.worker , args = (i,inq, retq )  ) }
			#processes[i]['process'].daemon = True
			processes[i]['process'].start()

		for batch in range(0, len(traindf) , chunksize ):

			slicedf = traindf.iloc[batch:batch+chunksize, :]
			fams = list(set(list(slicedf.HogFamA.unique()) + list(slicedf.HogFamB.unique() ) ) )
			total= {}

			for fam in fams:
				orthxml = self.READ_ORTHO(fam)
				if orthxml is not None:
					inq.put((fam,orthxml))
			done = []
			count = 0
			while len(fams)-1 > count:
				try:
					data =retq.get(False)
					count+=1
					total.update(data)
				except :
					pass
				time.sleep(.01)

			gc.collect()
			retdf= pd.DataFrame.from_dict( total , orient= 'index')
			slicedf = slicedf.merge( retdf , left_on = 'HogFamA' , right_index = True , how= 'left')
			slicedf = slicedf.merge( retdf , left_on = 'HogFamB' , right_index = True , how= 'left')
			slicedf = slicedf.dropna(subset=['mat_y', 'mat_x'] , how = 'any')
			slicedf['xtrain'] = slicedf.apply( calculate_x , axis = 1)
			X_train = np.vstack( slicedf['xtrain'])
			y_train = slicedf.truth
			print(slicedf)

			yield (X_train, y_train)
		for i in processes:
			inq.put(None)
		for i in processes:
			processes[i]['process'].terminate()

	def retmat_mp_profiles(self, fams , nworkers = 25, chunksize=50 , verbose = False ):
		"""
		function used to create dataframe containing binary profiles
		and trees of fams
		"""

		fams = [ f for f in fams if f]
		retq= mp.Queue(-1)
		inq= mp.Queue(-1)
		processes = {}
		mp.log_to_stderr()
		logger = mp.get_logger()
		logger.setLevel(logging.INFO)
		total = {}

		for i in range(nworkers):
			processes[i] = {'time':time.time() , 'process': mp.Process( target = self.worker , args = (i,inq, retq )  ) }
			#processes[i]['process'].daemon = True
			processes[i]['process'].start()
		for fam in fams:
			if verbose == True:
				print(fam)
			try:
				orthxml = self.READ_ORTHO(fam)
			except:
				orthxml = None
			if orthxml is not None:
				inq.put((fam,orthxml))
		done = []
		count = 0

		while len(fams)-1 > count :
			try:
				data =retq.get(False	)
				count+=1
				total.update(data)
				if count % 100 == 0 :
					print(count)
			except :
				pass
			time.sleep(.01)

		for i in range(nworkers):
			processes[i]['process'].terminate()
		retdf= pd.DataFrame.from_dict( total , orient= 'index')
		return retdf

	def hog_query(self, hog_id=None, fam_id=None , k = 100  ):
		"""
		Given a hog_id or a fam_id as a query, returns a dictionary containing the results of the LSH.
		:param hog_id: query hog id
		:param fam_id: query fam id
		:return: list containing the results of the LSH for the given query
		"""

		if hog_id is not None:
			fam_id = self.hogid2fam(hog_id)
		query_hash = hashutils.fam2hash_hdf5(fam_id, self.hashes_h5 , nsamples=  self.nsamples )
		#print(query_hash.hashvalues)
		results = self.lshobj.query(query_hash, k)


		return results

	def hog_query_sorted(self, hog_id=None, fam_id=None , k = 100  ):
		"""
		Given a hog_id or a fam_id as a query, returns a dictionary containing the results of the LSH.
		:param hog_id: query hog id
		:param fam_id: query fam id
		:return: list containing the results of the LSH for the given query
		"""

		if hog_id is not None:
			fam_id = self.hogid2fam(hog_id)
		query_hash = hashutils.fam2hash_hdf5(fam_id, self.hashes_h5 , nsamples=  self.nsamples )
		results = self.lshobj.query(query_hash, k)
		hogdict = self.pull_hashes(results)

		hogdict = { hog: hogdict[hog].jaccard(query_hash) for hog in hogdict  }
		sortedhogs = [(k, v) for k, v in hogdict.items()]
		sortedhogs = sorted(student_tuples, key=lambda x: x[1])
		sortedhogs = [ h[0] for h in sortehogs.reverse() ]
		return hogdict

	def pull_hashes(self , hoglist):

		"""
		Given a list of hog_ids , returns a dictionary containing their hashes.
		This uses the hdf5 file to get the hashvalues
		:param hog_id: query hog id
		:param fam_id: query fam id
		:return: a dict containing the hash values of the hogs in hoglist
		"""

		return { entry: hashutils.fam2hash_hdf5( self.hogid2fam(entry), self.hashes_h5 , nsamples=  self.nsamples) for entry in hoglist}

	def pull_matrows(self,fams):
		"""
		given a list of fams return the submatrix containing their profiles

		:return:fams sorted, sparse mat
		"""
		return self.profile_matrix[np.asarray(fams),:]


	@staticmethod
	def sort_hashes(query_hash,hashes):
		"""
		Given a dict of hogs:hashes, returns a sorted array of hogs and jaccard distances relative to query hog.
		:param query hash: weighted minhash of the query
		:param hashes: a dict of hogs:hashes
		:return: sortedhogs, jaccard
		"""
		#sort the hashes by their jaccard relative to query hash
		jaccard=[ query_hash.jaccard(hashes[hog]) for hog in hashes]
		index = np.argsort(jaccard)
		sortedhogs = np.asarry(list(hashes.keys()))[index]
		jaccard= jaccard[index]
		return sortedhogs, jaccard

	@staticmethod
	def allvall_hashes(hashes):
		"""
		Given a dict of hogs:hashes, returns generate an all v all jaccard distance matrix.
		:param hashes: a dict of hogs:hashes
		:return: hashmat
		"""
		#generate an all v all jaccard distance matrix
		hashmat = np.zeros((len(hashes),len(hashes)))
		for i , hog1 in enumerate(hashes):
			for j, hog2 in enumerate(hashes):
				if i < j :
					hashmat[i,j]= hashes[hog1].jaccard(hashes[hog2])
		hashmat = hashmat+hashmat.T
		np.fill_diagonal(hashmat, 1)
		return hashmat

	def hog_v_hog(self, hogs):
		"""
		give two hogs returns jaccard distance.
		:param hog1 , hog2: str hog id
		:return: jaccard score
		"""
		hog1,hog2 = hogs
		#generate an all v all jaccard distance matrix
		hashes = self.pull_hashes([hog1,hog2])
		hashes = list(hashes.values())
		return hashes[0].jaccard(hashes[1])

	def allvall_nx(G,hashes,thresh =None):

		"""
		Given a dict of hogs:hashes, returns generate an all v all jaccard distance matrix.
		:param hashes: a dict of hogs:hashes
		:return: hashmat
		"""

		#generate an all v all jaccard distance matrix

		hashmat = [[ hashes[hog1].jaccard(hashes[hog2]) if j>i else 0 for j,hog2 in enumerate(hashes[0:i] ) ] for i,hog1 in enumerate(hashes) ]
		hashmat = np.asarray(hashmat)
		hashmat+= hashmat.T
		np.fill_diagonal(hashmat, 1)

		#hashmat = np.zeros((len(hashes),len(hashes)))

		#for i , hog1 in enumerate(hashes):
		#	for j, hog2 in enumerate(hashes):
		#		hashmat[i,j]= hashes[hog1].jaccard(hashes[hog2])
		return hashmat

	def iternetwork(seedHOG):
		pass

	def rank_hashes(query_hash,hashes):
		jaccard = []
		sorted = []
		scores = {}
		hogsRanked = np.asarray(list(hashes.keys()))
		for i, hog in enumerate(hashes):
			score = query_hash.jaccard(hashes[hog])
			jaccard.append( score)
			scores[hog] = score
		hogsRanked = list( hogsRanked[ np.argsort(jaccard) ] )
		jaccard = np.sort(jaccard)
		return hogsRanked, jaccard
