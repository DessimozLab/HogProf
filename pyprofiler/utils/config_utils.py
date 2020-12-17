import json
import os

#don't forget the trailing slash on all dirs

config = {
    "dir":{
#    "datadir": "/home/cactuskid13/mntpt/unil_backup/profilingbackup/",
    "datadir": "/home/cactuskid13/pyprofiler/pyprofiler/plants/",

    "omadir": "/home/cactuskid13/mntpt/OMA/latest/"
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
