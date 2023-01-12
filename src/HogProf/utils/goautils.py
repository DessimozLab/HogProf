
from __future__ import print_function

from goatools import semantic
from pyoma.browser.models import ProteinEntry
from goatools.obo_parser import GODag
from . import hashutils
import pickle
from goatools.go_enrichment import GOEnrichmentStudy
from pyoma.browser import db
import multiprocessing as mp
from tables import *
import time
##############enrichment##############################################

def return_enrichment_study_obj(gaf_taxfiltered, obo = None):
    '''
    Generate go enrichment study object with a background dataset.
    '''
    if obo is None:
        obodag = GODag(config_utils.datadir+"/GOData/go-basic.obo")
    else:
        obodag = GODag(obo)

    goeaobj = GOEnrichmentStudy(
        gaf_taxfiltered.keys(), #
        gaf_taxfiltered, # geneid/GO associations possible with tree used for DB
        obodag, # Ontologies
        propagate_counts = False,
        alpha = 0.15, # default significance cut-off
        methods = ['fdr_bh']) # defult multipletest correction method
    return goeaobj


def buildGAF(gaf_file , universe= None):

    gaf_filtered = {}
    with open(gaf_file, mode='r') as gafin:
        for line in gafin:
            words = line.split()
            if words[0] not in gaf_filtered:
                gaf_filtered[words[0]]=set([words[1]])
            else:
                gaf_filtered[words[0]].add(words[1])

    if universe:
        gaf_filtered = { prot:gaf_filtered[prot] for prot in universe}


    return gaf_filtered


def return_hogs_timeout( result, retq):
    print('started')
    with open_file(config_utils.config['dir']['omadir']+'OmaServer.h5' , mode="r") as h5_oma:
        db_obj = db.Database(h5_oma)
        res =  [ ProteinEntry(db_obj, e).omaid for  e in db_obj.member_of_fam(int(result)) ]
        retq.put(res)

def run_GOEA_onresults(results, db_obj, goeaobj, outname = None):
    '''
        Perform enrichment analysis on returned results
        grabs all member protein of all hogs in result
        returns goe results and HOG composition
    '''
    hogids =[ "HOG:" + (7-len(str(fam_id))) * '0' + str(fam_id) for fam_id in results ]
    HOGS={}
    print('compiling hogs')
    prots = []

    print('mod13')
    retq = mp.Queue()

    for i,result in enumerate(results):
        if i %10 ==0:
            print(i)
        print(result)
        HOGS[result]=[]
        p = mp.Process( target= return_hogs_timeout , args= (result, retq))
        p.start()

        t0 = time.time()
        timeout = False

        while time.time()-t0 < 10 :
            time.sleep(.1)
            if p.is_alive() == False:
                print('done')
                break
        if time.time()-t0 > 10:
            timeout = True
            print('Dead')
        p.terminate()
        del p


        if retq.empty() == False:
            iterobj = retq.get(10)
            #retq get
            for k,member in enumerate(iterobj):
                if k < 1:
                    print(member)
                if k > 500:
                    break
                HOGS[result].append(member)
                prots.append(member)

    print('done')
    print('running GO enrichment study')
    goea_results_all = goeaobj.run_study(prots )
    print('done')
    if outname:
        with open( config_utils.datadir + outname + 'Hogs2Prots.pkl' , 'wb' ) as save:
            save.write(pickle.dumps(HOGS,2))
        goeaobj.wr_txt(config_utils.datadir+ str(outname)+"enrichment.txt", goea_results_all)
    print('DONE!')
    return goea_results_all, HOGS


def run_GOEA_onresults_tar(results, tar, goeaobj, outname = None):
    '''
        Perform enrichment analysis on returned results
        grabs all member protein of all hogs in result
        returns goe results and HOG composition
    '''
    ## TODO: finish this function with tar hog to list of prot IDS
    #print(db_obj.member_of_hog_id(int(results[0])))
    #hogids =[ "HOG:" + (7-len(fam_id)) * '0' + fam_id for fam_id in results ]
    #print( db_obj.member_of_hog_id(hogids[0]) )


    HOGS={}
    print('compiling hogs')
    prots = []
    for i,result in enumerate(hogids):
        if i %10 ==0:
            print(i)
        HOGS[result]=[]
        for member in db_obj.iter_members_of_hog_id(result):
            HOGS[result].append(member.omaid)
            prots.append(member.omaid)
    print('done')
    print('running GO enrichment study')

    goea_results_all = goeaobj.run_study(prots )
    print('done')
    with open( config_utils.datadir + outname + 'Hogs2Prots.pkl' , 'wb' ) as save:
       save.write(pickle.dumps(HOGS,2))

    goeaobj.wr_txt(config_utils.datadir+ str(outname)+"enrichment.txt", goea_results_all)
    print('DONE!')
    return goea_results_all, HOGS


######################resnik semsim ###################################################

def resnik_sim_hdf5(go_id1, go_id2, godag, termcounts, hdf5):
    '''
        Computes Resnik's similarity measure.
    '''
    try:
        msca_goid = deepest_common_ancestor_hdf5([goterm2id(go_id1), goterm2id(go_id2)], godag, hdf5)
        score = semantic.get_info_content(msca_goid, termcounts)
    except:
        score = -1
    return score


def deepest_common_ancestor_hdf5(go_ids, godag, hdf5):
    '''
        Gets the nearest common ancestor
        using the above function.
        Only returns single most specific - assumes unique exists.
    '''
    # Take the element at maximum depth.
    return max(common_parent_go_ids_hdf5(go_ids, hdf5), key=lambda t: godag[t].depth)

def common_parent_go_ids_hdf5(go_ids, hdf5_set):
    '''
        Finds the common ancestors in the GO
        tree of the list of goids in the input.
    '''
    candidates = set(hdf5_set[go_ids[0]].tolist())
    for go_id in go_ids[1:]:
        candidates_to_add = set(hdf5_set[go_id].tolist())
        candidates.intersection_update(candidates_to_add)
    corrected_candidates = [id2goterm(c) for c in candidates]
    return corrected_candidates

def resnik_sim_pandas(tup, df , termcounts):
    '''
        Computes Resnik's similarity measure.
    '''
    go_id1, go_id2 = tup
    #print(df.head())
    if go_id1 == go_id2:
        return semantic.get_info_content(go_id1, termcounts)

    elif go_id2 in df.index and go_id1 in df.index:

        ancestors = df.loc[str(go_id2)].parents
        ancestors += df.loc[str(go_id1)].parents
        terms = df.loc[ancestors]
        ancestors_set = terms.parents.tolist()
        intersection = set(ancestors_set[0]).intersection(* ancestors_set[1:])
        common_ancestors = df.loc[list(intersection)]
        common_ancestors = common_ancestors.sort_values('depth', ascending= False)
        msca_goid = common_ancestors.index.tolist()[0]
        return semantic.get_info_content(msca_goid, termcounts)

    else:
        return -1


def get_go_terms_gaf(hog_id, pyoma_dbobj, gaf , genomes = None):
    '''
        iterate over hog members and get the go information from a gaf in memory
    '''
    fam = hashutils.hogid2fam(hog_id)
    go_terms = { mr.omaid:gaf[mr.omaid] for mr in pyoma_dbobj.iter_members_of_hog_id(hog_id) if mr.omaid in gaf  }
    return go_terms


def goterm2id(go_term_to_modif):

    return int(go_term_to_modif.split(':')[1])

def id2goterm(go_term_to_modif):
    return 'GO:{:07d}'.format(go_term_to_modif)