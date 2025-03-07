from pyoma.browser import db
import pickle
import pandas as pd
import h5py
import random
from tables import *
import numpy as np
import random
import ete3
import argparse
import sys
#from validation import validation_semantic_similarity
from HogProf.utils import hashutils , pyhamutils , files_utils
from time import time
import multiprocessing as mp
import functools
import numpy as np
import time
import gc
import logging
import os
from pyoma.browser import db
np.random.seed(0)
random.seed(0)
class Profiler:

	"""
	A profiler object allows the user to query the LSH with HOGs and get a list of result HOGs back

	"""
	def __init__(self,lshforestpath = None, hashes_h5=None, mat_path= None, oma = False , nsamples = 256 , 
			  mastertree = None , reformat_names = False , swap2taxcode = False , use_phyloxml = False , 
			  taxfilter = None , taxmask = None, slicesubhogs = False ):
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
		slicesubhogs (bool): Whether to slice subhogs.

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
		self.slicesubhogs = slicesubhogs
		hamfunction = pyhamutils.get_ham_treemap_from_row
		hashfunction = hashutils.row2hash
		### special operations for sliced subhogs
		if self.slicesubhogs is True:
			self.fam2orthoxmlpath = os.path.join(os.path.dirname(lshforestpath), 'fam2orthoxml.csv')
			### adjust functions
			hamfunction = pyhamutils.get_subhog_ham_treemaps_from_row
			hashfunction = hashutils.hash_trees_subhogs
			### get a dictionary of subhog ids
			id2famsubhog_df = pd.read_csv(self.fam2orthoxmlpath) ### Athina note: here used to be index_col=0 !!!!!!!!!!!!!!!!!!!!!!!!
			# print for fam=15
			print(id2famsubhog_df.head())
			print(id2famsubhog_df[id2famsubhog_df['fam'] == 15]) 
			# Group by 'fam' and create a dictionary of indices
			fam_dict = id2famsubhog_df.groupby('fam', group_keys=False).apply(lambda x: x.index.tolist(), include_groups=False).to_dict()
			self.fam_dict = fam_dict
			# Create the reverse dictionary too
			subhogid_to_fam_dict = {subhogid: fam for fam, subhogid_list in fam_dict.items() for subhogid in subhogid_list}
			self.subhogid_to_fam_dict = subhogid_to_fam_dict
			# Reconstruct the subhog_id for each row and create a subhog_dict
			subhog_dict = id2famsubhog_df.apply(lambda x: f"{x['fam']}_{x['subhog_id']}", axis=1).to_dict()
			self.subhog_dict = subhog_dict
			# Create the reverse dictionary too
			subhogname_to_id_dict = {subhogname: subhogid for subhogid, subhogname in subhog_dict.items()}
			self.subhogname_to_id_dict = subhogname_to_id_dict
			# Connect subhognames to fams
			subhogname_to_fam_dict = {subhog_id: subhog_id.split('_')[0] for subhog_id in subhogname_to_id_dict.keys()}
			self.subhogname_to_fam_dict = subhogname_to_fam_dict
			# Connect subhogids to levels
			subhog_to_level_dict = {subhog_id: subhog_id.split('_')[1] for subhog_id in subhogname_to_id_dict.keys()}
			self.subhog_to_level_dict = subhog_to_level_dict
			# Connect HOG IDs to fams
			hogid_to_fam_dict = {subhogname.split("_")[-2]: fam for subhogname, fam in subhogname_to_fam_dict.items()}
			self.hogid_to_fam_dict = hogid_to_fam_dict
			
		print('h5' , self.hashes_h5 , self.hashes_h5.keys())
		self.nsamples = nsamples
		if 'xml' in mastertree.lower():
			project = Phyloxml()
			project.build_from_file(mastertree)
			trees = [t for t in  project.get_phylogeny()]
			self.tree = [ n for n in trees[0] ][0]
			self.use_phyloxml = True
			print('using phyloxml')
			print( 'loaded tree:' , self.tree )
			self.tree_string = mastertree
		else:
			try:
				self.tree = ete3.Tree(mastertree, format=1 , quoted_node_names= True)
				print( 'loaded tree:', self.tree )
			except:
				self.tree = ete3.Tree(mastertree, format=0)
		with open(mastertree) as treein:
			self.tree_string = treein.read()
		#self.tree_string = self.tree_ete3.write(format=0)

        
		if oma:
			self.reformat_names = reformat_names
			self.saving_path = mat_path
			if not os.path.exists(self.saving_path):
				os.makedirs(self.saving_path)
			
			if self.reformat_names:
				self.tree, self.idmapper = pyhamutils.tree2numerical(self.tree)
				self.tree_string = self.tree.write(format=1)
				with open( self.saving_path + 'reformatted_tree.nwk', 'w') as treeout:
					treeout.write(self.tree.write(format=0 ))
				with open( self.saving_path + 'idmapper.pkl', 'wb') as idout:
					idout.write( pickle.dumps(self.idmapper))
				print('reformatted tree')
				print( self.tree )
				self.tree_string = self.tree.write(format=1) 
				#remap taxfilter and taxmask
				if taxfilter:
					self.tax_filter = [ self.idmapper[tax] for tax in taxfilter ]
				if taxmask:
					self.tax_mask = self.idmapper[taxmask]
			else:
				self.idmapper = None
				self.tax_filter = taxfilter
				self.tax_mask = taxmask
				self.tree_string = self.tree.write(format=1)
			
			self.taxaIndex, self.ReverseTaxaIndex = files_utils.generate_taxa_index(self.tree)
			self.treeweights = hashutils.generate_treeweights(self.tree , self.taxaIndex , None, None )
			self.swap2taxcode = swap2taxcode
			self.use_phyloxml = use_phyloxml
			self.tax_filter = None
			self.tax_mask = None
			

			h5_oma = open_file(oma, mode="r")
			self.db_obj = db.Database(h5_oma)
			#self.treeweights = hashutils.generate_treeweights(self.tree , self.taxaIndex , None, None )
			

			self.READ_ORTHO = functools.partial(pyhamutils.get_orthoxml_oma, db_obj=self.db_obj)
			
			self.HAM_PIPELINE = functools.partial( hamfunction, tree=self.tree_string ,  swap_ids=self.swap2taxcode , reformat_names = self.reformat_names , 
												  orthoXML_as_string = True , use_phyloxml = self.use_phyloxml , orthomapper = self.idmapper ) 
			
			self.HASH_PIPELINE = functools.partial(hashfunction , taxaIndex=self.taxaIndex  , treeweights=self.treeweights , wmg=None )
		
		self.taxaIndex, self.ReverseTaxaIndex = files_utils.generate_taxa_index(self.tree)

		print('DONE with profiler init')

	def hogid2fam(self, hog_entry, reverse = False):
		#print("hog_entry", hog_entry)
		if type(hog_entry )== int and not self.slicesubhogs:
			return hog_entry
		elif self.slicesubhogs:
			### if there is a column without label, use it as index
			id2famsubhog_df = pd.read_csv(self.fam2orthoxmlpath, index_col=0)
			### check if index is labeled
			if id2famsubhog_df.index.name is not None:
				id2famsubhog_df = pd.read_csv(self.fam2orthoxmlpath)
			#print(id2famsubhog_df)
			fam_dict = id2famsubhog_df.groupby('fam', group_keys=False).apply(lambda x: x.index.tolist(), include_groups=False).to_dict()
			#subhog_dict = id2famsubhog_df.set_index('subhog_id').to_dict(orient='index')
			if isinstance(hog_entry, int):
				indices = fam_dict[hog_entry]
				if reverse:
					### return also the indices to subhogids dict
					#print(id2famsubhog_df)
					#subhog_dict = id2famsubhog_df.set_index('subhog_id').to_dict(orient='index')
					### for each row of the df with index in indices, reconstruct the subhog id from columns fam and subhog_id
					subhog_dict = id2famsubhog_df.loc[indices, ['fam', 'subhog_id']].apply(lambda x: f"{x['fam']}_{x['subhog_id']}", axis=1).to_dict()
					#print(subhog_dict)
					return indices, subhog_dict
			elif isinstance(hog_entry, str):
				fam_parts = hog_entry.split('_')
				fam_int = int(fam_parts[0])
				subhog_str = '_'.join(fam_parts[1:])
				indices = id2famsubhog_df[(id2famsubhog_df['fam'] == fam_int) & (id2famsubhog_df['subhog_id'] == subhog_str)].index.tolist()
				return indices, {indices[0]:hog_entry}
			else:
				return indices
		else:
			hog_entry = self.db_obj.entry_by_entry_nr(self.db_obj.id_resolver.resolve(hog_entry))
			famnr = int(self.db_obj.hog_family( entry=hog_entry ) )
			return famnr

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
				# Collect indices for losses, duplications, and presence of genes
				losses = [ self.taxaIndex[n.name]  for n in tp.traverse() if n.lost and n.name in self.taxaIndex  ]
				dupl = [ self.taxaIndex[n.name]  for n in tp.traverse() if n.dupl  and n.name in self.taxaIndex  ]
				presence = [ self.taxaIndex[n.name]  for n in tp.traverse() if n.nbr_genes > 0  and n.name in self.taxaIndex  ]
				# Create a dictionary mapping event types to their indices
				indices = dict(zip (['presence', 'loss', 'dup'],[presence,losses,dupl] ) )
				# Initialize a zero matrix
				hog_matrix_raw = np.zeros((1, 3*len(self.taxaIndex)))
				# Iterate over the indices dictionary and update the matrix if there are any indices for an event type
				for i,event in enumerate(indices):
					if len(indices[event])>0:
						taxindex = np.asarray(indices[event])
						hogindex = np.asarray(indices[event])+i*len(self.taxaIndex)
						hog_matrix_raw[:,hogindex] = 1
				# Put the result in the return queue
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

	def hog_query(self, hog_id=None, fam_id=None , k = 100):
		"""
		Given a hog_id or a fam_id as a query, returns a dictionary containing the results of the LSH.
		:param hog_id: query hog id
		:param fam_id: query fam id
		:return: list containing the results of the LSH for the given query
		"""

		if hog_id is not None:
			if not self.slicesubhogs:
				fam_id = hog_id
		if self.slicesubhogs:
			query_hash = hashutils.fam2hash_hdf5(fam_id, self.hashes_h5 , nsamples=  self.nsamples, fam2orthoxmlpath = self.fam2orthoxmlpath )
		else:
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
			if  self.slicesubhogs is False:
				fam_id = self.hogid2fam(hog_id)
			else:
				### if a family is query
				if isinstance(hog_id, int):
					fam_id = hog_id
					family_subhogids_list = self.fam_dict[hog_id]
					family_subhogs_list = [self.subhog_dict[i] for i in family_subhogids_list]
				### if a family is query by HOG ID
				elif 'HOG' in hog_id:
					fam_id = self.hogid_to_fam_dict[hog_id]
					family_subhogids_list = self.fam_dict[int(fam_id)]
					family_subhogs_list = [self.subhog_dict[i] for i in family_subhogids_list]
				### if a subhog is query
				else:
					fam_id = int(hog_id.split('_')[0])
					family_subhogs_list = [hog_id]
					family_subhogids_list = [self.subhogname_to_id_dict[hog_id]]

				#fam_id, subhog_dict = self.hogid2fam(hog_id, reverse=True)
		#if self.slicesubhogs and len(fam_id) == 1:
		#		fam_id = fam_id[0]
		### case where a family is given as query (including many subfamilies)
		if self.slicesubhogs:
			print("\nHOG query sorted - Case of sliced subhogs")
			#print(hog_id,'first fam_id', fam_id[0]) # family, [index1, index2, ...]
			query_hashes_dict = {i:hashutils.fam2hash_hdf5(i, self.hashes_h5 , nsamples=  self.nsamples ) for i in family_subhogids_list}
			#print('query_hashes',list(query_hashes_dict.keys())[0],query_hashes_dict[list(query_hashes_dict.keys())[0]]) # index1: hash1, index2: hash2, ...
			### works without filtering
			results_dict = { fam: self.lshobj.query(query_hash, k) for fam, query_hash in query_hashes_dict.items() }
			#print('results',list(results_dict.keys())[0],'first result hogid',results_dict[list(results_dict.keys())[0]][0]) # index1: results1, index2: results2, ...
			### the line below is not optimised!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			#print(query_hashes_dict)
			
			hogdict = { fam: self.pull_hashes(results) for fam, results in results_dict.items() }
			#print('hogdict',list(hogdict.keys())[0],'first result: 0_0_0',hogdict[list(hogdict.keys())[0]]['0_0_0']) # {index1:{hitsubhog1:hash1, hitsubhog2:hash2, ...}, index2:{subhog1:hash1, subhog2:hash2, ...}, ...}				
			hogdict = { fam: { subhog: hogdict[fam][subhog].jaccard(query_hashes_dict[fam]) for subhog in hogdict[fam] } for fam in hogdict }
			#print('hogdict',list(hogdict.keys())[0],'first result: 0_0_0',hogdict[list(hogdict.keys())[0]]['0_0_0']) # {index1: jaccard1, index2: jaccard2, ...}

			### making df instead of sorted list to include family info
			data = []
			for index, hit_dict in hogdict.items():
				for subhoghit, jaccard in hit_dict.items():
					data.append({
						'query_subhogid': self.subhog_dict[index],
						'hit_subhogid': subhoghit,
						'jaccard': jaccard
					})
			subhogsdf = pd.DataFrame(data, columns=['query_subhogid', 'hit_subhogid', 'jaccard'])
			subhogsdf['query_family'] = fam_id
			#print(subhogsdf.head())
			# Add query_family and query_level columns using the dictionaries
			subhogsdf['query_level'] = subhogsdf['query_subhogid'].map(self.subhog_to_level_dict)
			# Add hit_family and hit_level columns using the dictionaries
			subhogsdf['hit_family'] = subhogsdf['hit_subhogid'].map(self.subhogname_to_fam_dict)
			subhogsdf['hit_level'] = subhogsdf['hit_subhogid'].map(self.subhog_to_level_dict)
			subhogsdf = subhogsdf.sort_values(['jaccard', 'query_subhogid', 'hit_family', 'hit_level'], ascending=False)
			return hogdict, subhogsdf	
		else:
			query_hash = hashutils.fam2hash_hdf5(fam_id, self.hashes_h5 , nsamples=  self.nsamples)
			results = self.lshobj.query(query_hash, k)
			hogdict = self.pull_hashes(results)
			hogdict = { hog: hogdict[hog].jaccard(query_hash) for hog in hogdict  }
		#print(hogdict)
		sortedhogs = [(k, v) for k, v in hogdict.items()]
		sortedhogs = sorted(student_tuples, key=lambda x: x[1])
		sortedhogs = [ h[0] for h in sortehogs.reverse() ]
		return hogdict , sortedhogs

	def pull_hashes(self , hoglist):

		"""
		Given a list of hog_ids , returns a dictionary containing their hashes.
		This uses the hdf5 file to get the hashvalues
		:param hog_id: query hog id
		:param fam_id: query fam id
		:return: a dict containing the hash values of the hogs in hoglist
		"""
		#print('pulling hashes')
		### hogid2fam returns the index of the family in the h5 file for normal cases
		### for sliced subhogs, it returns the list of indices of the subhogs in the h5 file
		### returns indices, subhog_dict
		if not self.slicesubhogs:
			return { entry: hashutils.fam2hash_hdf5( self.hogid2fam(entry), self.hashes_h5 , nsamples=  self.nsamples) for entry in hoglist}
		else:
			result_dict = {}
			for entry in hoglist:
				indices, subhog_dict = self.hogid2fam(entry)
				result_dict[entry] = [hashutils.fam2hash_hdf5(idx, self.hashes_h5, nsamples=self.nsamples) for idx in indices][0]
				
			return result_dict

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
		pas

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

def main():
	parser = argparse.ArgumentParser()
	parser.add_argument('--lshforestpath', help='Path to the newlshforest.plk file',type = str)
	parser.add_argument('--hashes_h5', help='Path to the hashes.h5 file',type = str)
	parser.add_argument('--mat_path', help='Path to the fam2orthoxml.csv file' , type = str)
	parser.add_argument('--OMA', help='Path to OmaServer.h5 file' , type = str)
	parser.add_argument('--nsamples', help='no of samples' , type = int)
	parser.add_argument('--slicesubhogs', help='Use slicesubhogs mode' , type = str)
	parser.add_argument('--mastertree', help='master taxonomic tree. nodes should correspond to orthoxml' , type = str)
	args = vars(parser.parse_args(sys.argv[1:]))
	print("\nReading arguments")

	if args['slicesubhogs'] == 'True':
		args['slicesubhogs'] = True
	elif args['slicesubhogs'] == 'False':
		args['slicesubhogs'] = False

	p = Profiler(lshforestpath=args['lshforestpath'],
             hashes_h5=args['hashes_h5'],
             mat_path=args['mat_path'],
             oma=args['OMA'],
             nsamples=args['nsamples'],
			 slicesubhogs=args['slicesubhogs'],
             mastertree=args['mastertree']
			 )
	print("\nProfiler object created")
	#hogdict, sortedhogs = p.hog_query_sorted( hog_id= '0_0_0' , k = 20 )
	#hogdict, sortedhogs = p.hog_query_sorted( hog_id= 0 , k = 20 )
	hogdict, sortedhogs = p.hog_query_sorted( hog_id= 712231 , k = 20 )
	print(sortedhogs)
	print()
	nonself_samelevel_hits = filtered_hogs = sortedhogs[(sortedhogs['query_level'] == sortedhogs['hit_level']) &
													 (sortedhogs['query_subhogid'] != sortedhogs['hit_subhogid'])]
	print(nonself_samelevel_hits)
	outfile = os.path.join(os.path.dirname(args['lshforestpath']), 'nonself_samelevel_hits.csv')
	nonself_samelevel_hits.to_csv(outfile, index=False)
	print()
	if sortedhogs['hit_subhogid'].str.contains('3584_HOG:E0712231_0').any(): ## for local was 0_0_0, for curnagl was 0_4_0
		print('got hit!\n')
	else:
		print('Warning! Did not find itself!\n')


if __name__ == '__main__':
    main()