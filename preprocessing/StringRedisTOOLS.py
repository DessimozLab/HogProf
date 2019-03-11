import itertools
import redis

from utils import hashutils
from pyoma.browser import db
from time import time


def connect2IDmap():
    r1 = redis.StrictRedis(host='10.0.63.33', port=6379, db=0)
    return r1

def connect2Stringmap():
    r2 = redis.StrictRedis(host='10.0.63.33', port=6379, db=1)
    return r2

def clearDB(dbnum):
    r = redis.StrictRedis(host='10.0.63.33', port=6379, db=dbnum)
    r.flushdb()
    return r


def search(name, source_list):
    for x in source_list:
        if x['source'] == name:
            return x

def fam2stringID(dbobj, hog_id, r):
    members = dbobj.iter_members_of_hog_id(hog_id)
    oma_id_mapper = db.OmaIdMapper(dbobj)
    XrefIdMapper = db.XrefIdMapper(dbobj)
    xrefs = {oma_id_mapper.omaid_to_entry_nr(m.omaid):XrefIdMapper.map_entry_nr(oma_id_mapper.omaid_to_entry_nr(m.omaid)) for m in members}
    xrefs_source = {genome: search('SourceID', ids)['xref'] for genome, ids in xrefs.items()}
    print(xrefs_source)
    ret_string = {genome: string for genome, string in xrefs_source.items() if string is not None}
    return ret_string

def HOGvsHOG(allstring1, allstring2, r2, datapath):
    # protein1 protein2 neighborhood fusion cooccurence coexpression experimental database textmining combined_score
    refs = ['neighborhood', 'fusion', 'cooccurence', 'coexpression', 'experimental', 'database', 'textmining', 'combined_score']
    final = {}
    with open(datapath, 'r') as stringdata:
        for combos in itertools.product(allstring1, allstring2):
            id1, id2 = combos
            print('ids', id1, id2)
            # for id1, id2 in itertools.combinations(allstring1, allstring2):
            if r2.get(''.join(sorted([id1, id2]))) is not None:
                line = r2.get(''.join(sorted([id1, id2])))
                stringdata.seek(int(line), 0)
                words = stringdata.readline()
                # words = linecache.getline(datapath, int(line)).split()
                # words = stringdata.readline(int(line)).split()
                print(words)
                row_dict = {refs[i]: int(entry) for i, entry in enumerate(words[2:])}
                final[''.join(sorted([id1, id2]))] = row_dict
            else:
                print('nothing in db')
        return final

# map HOGS to OMA ID

# map OMIDS to UNIPROT

# MAP UNIPROT to string

# Find STRING vs STRING in files
