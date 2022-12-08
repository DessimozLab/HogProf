import json
import os

#don't forget the trailing slash on all dirs

config = {
    "dir":{
    "datadir": "/work/FAC/FBM/DBC/cdessim2/default/dmoi/datasets/birds/",
    "omadir": "/work/FAC/FBM/DBC/cdessim2/default/dmoi/datasets/birds/"
    },
    "orthoxmltar":"",
    "email": "dmoi@unil.ch"
}

datadir = config['dir']['datadir']
omadir = config['dir']['omadir']
email = config['email']
tarfile = config['orthoxmltar']

#ncbi ID in names ( weird orthoxml format for bird db)
ncbi_inID = False
mapdict_exceptions = {'ANAPL':'8839','CHICK':'9031','FICAL':'59894','JUNHY':'40217','MELGA':'9103','MELUD':'13146','PARMJ':'9157','SERCA':'9135' ,'TAEGU':'59729' }

if len(tarfile)==0:
    tarfile = None

for dir in config['dir'].values():
    if not os.path.isdir(dir):
        os.mkdir(path=dir)
print(config)
