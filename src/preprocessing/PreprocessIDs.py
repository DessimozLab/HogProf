from goatools import obo_parser
import numpy as np
import json
import pandas as pd

import sys
sys.path.insert(0, '..')
from utils import config_utils
from utils import goatools_utils
import pickle

"""
a script to turn a GO DAG into a dataframe object.
This is useful to perform quick resnik distance calculations
"""

if __name__ == '__main__':

    if preprocess_config.preprocessGO ==  True:
        #Preprocess all of the GO terms' parents to avoid looking at the DAG
        start_time = time.time()
        obo_reader = obo_parser.GODag(obo_file=config_utils.datadir + 'GOData/go-basic.obo')
        godict = {}
        for i,go_term in enumerate(obo_reader):
            go_term_read = obo_reader[go_term]
            if i %1000==0:
                print(go_term_read)
            if go_term_read.namespace == 'biological_process':
                godict[go_term_read.id]={}
                godict[go_term_read.id]['parents'] = [ p for p in go_term_read.get_all_parents()] + [go_term_read.id]
                godict[go_term_read.id]['level'] = go_term_read.level
                godict[go_term_read.id]['depth'] = go_term_read.depth
        goframe = pd.DataFrame.from_dict( godict , orient = 'index')
        print(goframe)
        open(config_utils.datadir + 'GOData/goframe.pkl' , 'wb').write(pickle.dumps( goframe, -1))
        print('Done with the parents in {} seconds'.format(time.time()-start_time))


print('DONE!')
