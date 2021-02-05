# Script for running IntaRNA on files which had no/not enough homologes
import subprocess as sp
import os
import time
import Bio.Seq as Seq

# path to the folder for IntaRNA
path = '/data/jannik/data/results/mafft/filtered/IntaRNA/'

# getting the names of the folders inside the IntaRNA directory
dirs = os.listdir(path)
start_time = time.perf_counter()
counter = 0

# looping over the directories
for dir in dirs:

	# change working directory to the folder containing the fasta file and
	# open it

	os.chdir(path+dir+'/')
	file = open(dir+'.fasta','r')

	# prepare the list for the sequence

	seq = []

	# transcribe the first entry (DNA of sRNA in c. glutamicum) to
	# corresponding RNA

	for entry in file:
		seq.append(entry)
	seq[1] = Seq.Seq(seq[1].replace('\n','')).complement().transcribe()
	file.close()

	# reopen file and write the RNA sequence plus identifier

	file = open(dir+'.fasta','w')
	file.write(seq[0])
	file.write(str(seq[1]))
	file.close()

	# run IntaRNA on the RNA and the mRNAs from the protein-coding genes; to
	# get the reference file, run transcribe_genes_to_mrna.py in the
	# scripts/sequences folder with your gff annotation file and genome

	command = 'IntaRNA -q '+ dir+'.fasta -t /data/jannik/data/sequences/refseq/mrna.fasta --threads 24 --outMode C  --outCsvCols id1,start1,end1,id2,start2,subseqDP,hybridDP,E ' + dir+'_intarna.csv --noSeed'
	process = sp.Popen(command.split(), stdout=sp.PIPE)
	output, error = process.communicate()
	print(output,'\n', error)
	counter += 1
	to_go = len(dirs)-counter
	time_used = (time.perf_counter()-start_time)/3600
	time_est = (time_used/counter)/24*to_go
	print('IntaRNA: Hours run: ' + str(time_used)+'est. Days to go: '+ str(time_est))
	print(str(counter)+' files ran, '+ str(len(dirs)-counter)+' to go!')
