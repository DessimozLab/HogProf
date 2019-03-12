import networkx as nx
import pandas as pd
from matplotlib import pyplot as plt
import glob
from pyoma.browser import db
import pickle



omadir = '/scratch/ul/projects/cdessimo/oma-browser/All.Jun2018/data'
db = db.Database(omadir + '/OmaServer.h5')


print('loading mapping')
experiments = ' fusion coexpression experiments textmining'
unidf = pd.read_csv( 'full_uniprot_2_string.04_2015.tsv' ,delim_whitespace = True , header = 0 )
unidf.columns = [ col.split('|')[0].replace('#','') for col in unidf.columns ]
unidf['uniprot_code'] = unidf.uniprot_ac.map(lambda x : x.split('|')[0])
unidf['uniprot_ac'] = unidf.uniprot_ac.map(lambda x : x.split('|')[1])
omadf = pd.read_csv( 'oma-uniprot.txt' ,delim_whitespace = True  , comment= '#' , names = ['oma' , 'uniprot'])
print('done')

print('loading network files')

networks = glob.glob('./*protein.links.full*txt')
print(networks)


def ID2HOG(ID):
    try:
        prot = db.id_resolver.resolve(ID)
        entry = db.entry_by_entry_nr(prot)
        HOG = entry[4]
    except:
        return None

    if len(HOG)>0:
        return HOG
    else:
        return None


for file in networks:
    ##species   uniprot_ac|uniprot_id   string_id   identity   bit_score
    #626523  C4G7Z8|C4G7Z8_9FIRM     GCWU000342_00164        100.00   608
    nxdf = pd.read_csv(file, sep=' ', header= 0)
    nxdf['prot1cut'] = nxdf.protein1.map( lambda x : x.split('.')[1])
    nxdf['prot2cut'] = nxdf.protein2.map( lambda x : x.split('.')[1])
    for col in nxdf.columns:
        if col+'_transferred' in nxdf.columns:
            nxdf[col]=nxdf[col].add( nxdf[col+'_transferred'] )
    nxdf['experimental'] = 0

    for data in experiments.split():
        nxdf['experimental']=nxdf['experimental'].add( nxdf[data.strip()])



    nxdf['uni1'] = nxdf.merge( unidf , how='left', left_on='prot1cut' , right_on='string_id' )['uniprot_code']
    nxdf['uni2'] = nxdf.merge( unidf , how='left', left_on='prot2cut' , right_on='string_id' )['uniprot_code']

    nxdf['oma1'] = nxdf.merge( omadf , how='left', left_on='uni1' , right_on='uniprot' )['oma']
    nxdf['oma2'] = nxdf.merge( omadf , how='left', left_on='uni2' , right_on='uniprot' )['oma']


    print(nxdf)
    count_nan = len(nxdf) - nxdf.count()
    print(count_nan)

    nxdfnotmapped = nxdf[ nxdf.oma1.isnull() | nxdf.oma2.isnull()]
    prots = set(list(nxdfnotmapped.uni1.unique())+list(nxdfnotmapped.uni2.unique()))
    hogmapUNI = { prot : ID2HOG(prot) for prot in prots}
    prots = set(list(nxdf.oma1.unique())+list(nxdf.oma2.unique()))
    hogmapOMA = { prot : ID2HOG(prot) for prot in prots}

    nxdf['HOG1'] = None
    nxdf['HOG2'] = None
    nxdf['HOG1'][nxdf.oma1.notnull()] =nxdf.oma1[nxdf.oma1.notnull()].map(hogmapOMA)
    nxdf['HOG2'][nxdf.oma2.notnull()] =nxdf.oma2[nxdf.oma2.notnull()].map(hogmapOMA)
    nxdf['HOG1'][nxdf.oma1.isnull()] =nxdf.uni1[nxdf.oma1.isnull()].map(hogmapUNI)
    nxdf['HOG2'][nxdf.oma2.isnull()] =nxdf.uni2[nxdf.oma2.isnull()].map(hogmapUNI)

    #leftovers
    prots = set(list(nxdf.prot1cut[nxdf.HOG1.isnull()].unique()) + list(nxdf.prot2cut[nxdf.HOG2.isnull()].unique() ))
    hogmapHAILMARRY = { prot : ID2HOG(prot) for prot in prots}
    nxdf['HOG1'][nxdf.HOG1.isnull()] =nxdf.prot1cut[nxdf.HOG1.isnull()].map(hogmapHAILMARRY)
    nxdf['HOG2'][nxdf.HOG2.isnull()] =nxdf.prot2cut[nxdf.HOG2.isnull()].map(hogmapHAILMARRY)

    count_nan = len(nxdf) - nxdf.count()
    print(count_nan)
    nxdf = nxdf[nxdf.HOG1.notnull() & nxdf.HOG2.notnull()]
    count_nan = len(nxdf) - nxdf.count()
    print(count_nan)

    nxdf.to_csv( file+'hogmap.csv')

    #read in protein links
    nxdf = nxdf[nxdf.experimental>300]
    print('builing graph')
    links = [ (row.HOG1, row.HOG2) for index,row in nxdf.iterrows()]
    #define graph
    G = nx.graph.Graph()
    G.add_edges_from(links)
    print(links[0:100])

    with open( file + 'links.pkl' , mode='wb') as graphout:
        graphout.write(pickle.dumps(links,-1))
    with open( file + '.pkl' , mode='wb') as graphout:
        graphout.write(pickle.dumps(G,-1))
    print('done graph')

print('DONE!')
