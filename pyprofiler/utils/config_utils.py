import json
import os


config = {
    "dir":{
    "datadir": "./profileDBs/",
    "omadir": "/scratch/cluster/monthly/dmoi/oma/june/"
    },
    "orthoxmltar":"",
    "email": "dmoi@unil.ch"
}

datadir = config['dir']['datadir']
omadir = config['dir']['omadir']
email = config['email']
tarfile = config['orthoxmltar']
if len(tarfile)==0:
    tarfile = None

for dir in config['dir'].values():
    if not os.path.isdir(dir):
        os.mkdir(path=dir)
print(config)
