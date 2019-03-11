import h5py
from goatools import obo_parser


def goterm2id(go_term_to_modif):
    id = int(go_term_to_modif.split(':')[1])
    return id


obo_iterator = obo_parser.OBOReader(obo_file='/home/laurent/Documents/project/phyloprofiling/data/go.obo')
obo_reader = obo_parser.GODag(obo_file='/home/laurent/Documents/project/phyloprofiling/data/go.obo')

with h5py.File('/home/laurent/Documents/project/phyloprofiling/data/parents.h5','r') as parents:
    numToTest = 456

    print(parents['go_terms'][numToTest])

    save = parents['go_terms'][numToTest].tolist()
    print(save)

    for go_term in obo_iterator.__iter__():
        if goterm2id(go_term.id) == numToTest:
            try:
                go_term_read = obo_reader[go_term.id]
                go_term_parents = go_term_read.get_all_parents()
                go_term_parents_int = [goterm2id(parent) for parent in go_term_parents]
                print(go_term_parents_int)
            except:
                print('BUG')
