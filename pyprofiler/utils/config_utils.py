import json
import os

with open( '../config.json' , 'r') as configin:
    config = json.loads(configin.read())
datadir = config['dir']['datadir']
omadir = config['dir']['omadir']
email = config['email']
for dir in config['dir'].values():
    if not os.path.isdir(dir):
        os.mkdir(path=dir)
print(config)
