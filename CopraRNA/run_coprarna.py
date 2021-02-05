# script for running CopraRNA for files

import os
import subprocess as sp
import time
start_time = time.perf_counter()
print('started:')
print(start_time)
path = '/data/jannik/data/results/mafft/filtered/CopraRNA/'
folders = os.listdir(path)
num_folders = len(folders)
counter = 0
for folder in folders:
	file_path = path + folder
	os.chdir(file_path)
	name = folder+'.fasta'
	print(name)
	command = '/home/jseidel/CopraRNA/CopraRNA-master/CopraRNA2.pl -srnaseq '+name+'  -websrv -topcount 200 -cores 36'
	process = sp.Popen(command.split(), stdout=sp.PIPE)
	output, error = process.communicate()
	print(output)
	print(error)
	diff_time = (time.perf_counter()-start_time)/3600
	counter += 1
	est_time = (diff_time/counter)/60*(num_folders-counter)
	print('CopraRNA:'+ str(diff_time)+'h since start of script, estimate ' +str(int(est_time)) +' d to go!')
	print(str(counter)+' of '+str(num_folders)+' finished! ('+str(counter/num_folders*100)+'%)')
	print('_________________________________________________________________________')
