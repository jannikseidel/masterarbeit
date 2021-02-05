# script for moving fasta-files in subfolders with the same name
# not needed if examin_filtered.py was run on the files
import os
import shutil

files = os.listdir('/data/jannik/data/results/mafft/filtered/IntaRNA/')

names = []
for file in files:
	if 'fasta' in file:
		name = file.replace('.fasta','')
		names.append(name)
for name in names:
	folder = '/data/jannik/data/results/mafft/filtered/IntaRNA/'+name
	print(folder)
	os.mkdir(folder)
	shutil.copy(folder+'.fasta',folder+'/'+name+'.fasta')
