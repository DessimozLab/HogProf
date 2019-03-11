from goatools import obo_parser
import h5py
import numpy as np
import tables
import json
from pyoma.browser import db

import pandas as pd
import time

import sys
sys.path.insert(0, '..')

from utils import config_utils
from utils import preprocess_config
from utils import goatools_utils
import pickle

def yield_hogs_with_annotations(annotation_dataset):
    for fam, annotations in enumerate(hogs):
        try:
            obj = json.loads(annotations)
            if obj and type(obj) is dict and len(obj) > 0:
                yield {fam : obj}
        except ValueError:
            pass

def goterm2id(go_term_to_modif):
    return_id = int(go_term_to_modif.split(':')[1])
    return return_id


def _get_entry_gene_ontology(oma, entry):
    return oma.root.Annotations.GeneOntology.read_where('EntryNr == {}'.format(entry))

def _get_hog_members(hog_id, oma):
    """
    Gets all gene members from the hog
    :param hog_id: hog id
    :return: list of genes from the hog
    """
    iterator = _iter_hog_members(hog_id, oma)
    population = frozenset([x['EntryNr'] for x in iterator])
    return population

def _hog_lex_range(hog):
    """
    Decodes hog format
    :param hog: (bytes or string): hog id
    :return: hog_str: encoded hog id
    """
    hog_str = hog.decode() if isinstance(hog, bytes) else hog
    return hog_str.encode()


def _iter_hog_members(hog_id, oma):
    """
    iterator over hog members / get genes
    :param hog_id: hog id
    :return: yields members of hog
    """
    hog_range = _hog_lex_range(hog_id)
    it = oma.root.Protein.Entries.where('({!r} <= OmaHOG) & (OmaHOG < {!r})'.format(*hog_range))
    for row in it:
        yield row.fetch_all_fields()

def iter_hog(oma):
    for hog in oma.root.HogLevel:
        yield hog[0]

def fam2hogid(fam_id):
    """
    Get hog id given fam
    :param fam_id:o fam
    :return: hog id
    """
    hog_id = "HOG:" + (7-len(str(fam_id))) * '0' + str(fam_id)
    return hog_id

def _get_go_terms(hog_id, oma, go_file):
    """
    Fetch the genes from hog id, then get all GO terms from it
    :param hog_id: hog id
    :return: list of GO term
    """
    genes = _get_hog_members(hog_id, oma)
    go_dict = {entry: {_format_go_term(e) for e in _get_entry_gene_ontology(oma, entry)} for entry in genes}
    go_dict_filtered = _filter_result(go_dict, go_file)

    return _clean_dictionary(go_dict_filtered)


def _format_go_term(e):
    # return e['TermNr']
    return 'GO:{:07d}'.format(e['TermNr'])


def _filter_result(go_dict, go_file):
    go_dict_filtered = {}
    for gene_name, terms in go_dict.items():
        filtered_terms = _filter_namespace(terms, go_file)
        if filtered_terms:
            go_dict_filtered[gene_name] = filtered_terms
    return go_dict_filtered

def _filter_namespace(list_terms, go_file, name='biological_process'):
    """
    Keep only go terms within the correct ontology
    :param list_terms: list of go terms
    :param name: namespace to keep; default: 'biological_process'
    :return: list of terms with the correct namespace
    """
    terms_to_keep = []
    for term in list_terms:
        try:
            if go_file[term].namespace == name:
                terms_to_keep.append(term)
        except KeyError:
            pass
    return terms_to_keep


def _clean_dictionary(dictionary):
    return {k: v for k, v in dictionary.items() if v}

if __name__ == '__main__':
    #if preprocess_config.preprocessGO ==True:
    if preprocess_config.preprocessGO ==  False:

        #Preprocess all of the GO terms' parents to avoid looking at the DAG
        obo_reader = obo_parser.GODag(obo_file=config_utils.datadir + 'GOData/go-basic.obo')
        dt = h5py.special_dtype(vlen=np.dtype('int32'))
        omah5 = tables.open_file(config_utils.omadir + 'OmaServer.h5', mode='r')
        with h5py.File(config_utils.datadir + 'GOData/GOparents.h5', 'w', libver='latest') as h5_go_terms:
            start_time = time()
            h5_go_terms.create_dataset('goterms2parents', (10000000,), dtype=dt)
            dataset_go_terms_parents = h5_go_terms['goterms2parents']
            count = 0
            for go_term in obo_reader:
                go_term_read = obo_reader[go_term]
                if go_term_read.namespace == 'biological_process':
                    go_term_parents = go_term_read.get_all_parents()

                    go_term_parents_int = [goterm2id(go_term_read.id)] + [goterm2id(parent) for parent in go_term_parents]
                    dataset_go_terms_parents[goterm2id(go_term_read.id)] = go_term_parents_int
                    count += 1
                    if count % 1000 == 0:
                        print('saving')
                        h5_go_terms.flush()
            h5_go_terms.flush()
            print('Done with the parents in {} seconds'.format(time()-start_time))

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
