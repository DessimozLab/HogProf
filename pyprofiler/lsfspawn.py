

import subprocess
import shlex
from utils import config_utils
import sys

#spawn lsf jobs to optimize db
#for loop with bo determined points
#try out the first few points with no slope


def openprocess(args , inputstr =None , verbose = False , wait = True):
	args = shlex.split(args)
	p = subprocess.Popen(args,  shell = True, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr= subprocess.PIPE )
	if verbose == True and inputstr is not None:
		print(inputstr.decode())

	if inputstr != None:
		p.stdin.write(inputstr.encode())
	if wait == True:
		output = p.communicate()
		if verbose == True:
			print(output)
		p.wait()
		return output[0].decode()
	else:
		return p

def bsub( script_path , test = False):
	argstr = 'bsub < '+script_path
	if test == True:
		print(argstr)
		script = open(script_path,'r').read()
		print(script)
	else:
		print(argstr)
		p = openprocess(argstr, inputstr =None, verbose = True , wait = True)
		return p



if __name__ == '__main__':

	lsfproto = open( './lsfproto.txt', 'r').read()
	names = [ 'pw' , 'lw', 'dw' ]
	params = [ [0,0,1], [0,1,0] , [1,0,0] ]
	dbtype = 'all'
	count = 0
	npoints = 10
	for i,p in enumerate(params):
		for n in names:
			 args = ' '.join([ ' --'+names[j]+ ' '+ str(p[j]) for j in range(len(p)) ])
		#create savedir for Run
		args += ' --dir '+config_utils.datadir + 'run_'+dbtype+str(count)+'/'
		args += ' --db ' +dbtype
		#create bash script
		script = lsfproto.replace( '[params]', args)
		script = script.replace( '[run]', str(count))
		#send it of with bsub
		new_script_path = './run_'+dbtype+str(count)+'.sh'
		with open( new_script_path , 'w') as scriptout:
			scriptout.write(script)
		count+=1

		bsub(new_script_path, test = False)


	for k in range(npoints):
		#just generate more npoints
		args = ' --dir '+config_utils.datadir + 'run_'+dbtype+str(count)+'/'
		args += ' --db ' +dbtype
		script = lsfproto.replace( '[params]', args)
		script = script.replace( '[run]', str(i))
		new_script_path = './run_'+dbtype+str(count)+'.sh'
		with open( new_script_path , 'w') as scriptout:
			scriptout.write(script)
		count+=1
		bsub(new_script_path, test = False)
