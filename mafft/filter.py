# script for filtering the mafft msa files for >=80% identity and assigning
# accession numbers to the sequences
# outputstatistics:
#			filename	positives	negatives	amount of Corynebacterium glutamicum

import os
import Bio
import Bio.Entrez as Entrez
import urllib
from matplotlib import pyplot as plt
import urllib3
import time
from fuzzywuzzy import fuzz
from fuzzywuzzy import process

# get new list of Organisms for CorpaRNA
url = 'http://rna.informatik.uni-freiburg.de/CopraRNA/CopraRNA_available_organisms.txt'
http = urllib3.PoolManager()
http_request = http.request('GET', url)

# decode data from CopraRNA list
data = http_request.data.decode('utf-8')
# splitting string at new line
data = data.split('\n')
counter = 0
new_data = []
# split each entry at tab and append to new list
for entry in data:
	if len(entry) < 1 or entry[0] == '#' or entry[0:3] == 'Ref':
		counter += 1
		continue
	entry = entry.split('\t')
	new_data.append(entry)
	counter += 1


ref_seq = {}
# reading the refseq accsessions with tax.name to dictionary
for entry in new_data:
	acc_num = []
	name =''
	for ent in entry:
		if ent[0] == 'N' and ent[1].isupper():
			for en in ent.split(' '):
				acc_num.append(en)
		else:
			name = ent.split(',')[0]
	if name not in ref_seq.keys():
		ref_seq[name] = acc_num
	else:
		for entry in ref_seq[name]:
			acc_num.append(entry)
		ref_seq[name] = acc_num

# specifing email for using NCBI Entrez
Entrez.email = 'st145304@stud.uni-stuttgart.de'
# get the file names of the fasta files from mafft and append to name list
files = os.listdir('/data/jannik/data/results/mafft/')
names=[]
for file in files:
	if '.fasta' not in file:
		continue
	if 'phylo' in file:
		continue
	else:
		names.append(file)
# check wether the script already ran for a certain ammount and remove analysed
# names from list, otherwise prepare a list for the analysed filenames
if os.path.isfile('/data/jannik/data/results/mafft/name_list.txt'):
	name_file = open('/data/jannik/data/results/mafft/name_list.txt','r')
	for entry in name_file:
		names.remove(entry.strip('\n'))
	name_file.close()
else:
	name_file = open('/data/jannik/data/results/mafft/name_list.txt','w')
	name_file.close()

# open the statistics file
stat= open('/data/jannik/data/results/mafft/filtered/stats.txt','w')

firstline = '#filename\tpositives\tnegatives\tc_glutamicum\n'
stat.write(firstline)
count_end = len(names)
stat_count = len(names)

# loop over the files from mafft output, while names try to do the the matching
# if error occures wait 60 sec and continue; used names are written to name_file
while names:
	try:
		for name in names:
			# open the file name list of proccesed files
			name_file = open('/data/jannik/data/results/mafft/name_list.txt','a')
			# path to the files
			path = '/data/jannik/data/results/mafft/'
			# combining path with file name
			path = path+name
			# open the file
			file = open(path,'r')
			# reading all the sequences into dictionary with the id as key
			id = ''
			sequences = {}
			for entry in file:
				if entry[0] == '>':
					sequences[entry[1:].strip('\n')]=''
					id = entry[1:].strip('\n')
				else:
					seq = sequences[id]
					seq = seq+entry.strip('\n')
					sequences[id] = seq
			# preparing the dictionaries for the classified sequences and
			# appending the sequences from C. glutamicum atcc 13032 to the
			# Dictionary
			pos_seq = {}
			pos_seq['NC_006958'] = sequences['NC_006958.1'].replace('-','')
			neg_seq = {}
			count_cor_glut = 0
			# loop over the sequences of the file
			for key in sequences.keys():
				ref = sequences['NC_006958.1']
				bound_start = 0
				bound_stop = len(ref)-1
				# if reference start or end begins with '-' search for first occurence
				# of a base; only this slice is used for calculating PI
				if ref[0] == '-':
					for i in range(0,len(ref)):
						bound_start = i
						if ref[i] != '-':
							break
				if ref[len(ref)-1] == '-':
					for i in reversed(range(0,len(ref))):
						bound_stop = i
						if ref[i] != '-':
							break

				# if the sequence is part of the reference sequence continue; else
				# calculate PI to ref-seq and if PI >= 80% search for the gi in the
				# nuccore NCBI Database, extract the taxonomic id, then search in the
				# taxonomy Database for the corresponding tax. name, if the name is
				# missing in the ref_seq keys, print this missing and continue,
				# else assign the first Accession id from ref_seq to a new variable
				# and append to the pos_seq dict with accession id as key

				if key == 'NC_006958.1' or key == 'NC_006958':
					continue
				else:
					counter = 0
					for i in range(bound_start,bound_stop):
						if ref[i] == sequences[key][i]:
							counter += 1
					if counter/(bound_stop-bound_start) >=0.8:
						new_key = key.split('|')
						new_key = [new_key[0],new_key[1],new_key[2],new_key[3]]
						final_key = ''
						summar = Entrez.esummary(db='nuccore',id = new_key[1])
						record = Entrez.read(summar)
						summar.close()
						summar1 = Entrez.esummary(db='taxonomy',id = record[0]['TaxId'])
						record1 = Entrez.read(summar1)
						summar1.close()
						tax_name = record1[0]['ScientificName']
						gi = new_key[1]
						db = [new_key[2],new_key[3]]
						if 'Corynebacterium glutamicum' in tax_name:
							count_cor_glut += 1
							continue
						# check if species is in ref_seq for CopraRNA;
						# use the string which matches the full tax name closest
						full_tax_name = tax_name
						tax_name = tax_name.split(' ')
						tax_name = tax_name[0] + '_' + tax_name[1]
						tax_names = []
						for key in ref_seq.keys():
							if key.startswith(tax_name):
								tax_names.append(key)
						print(tax_names)
						print(len(tax_names))
						if len(tax_names) == 1:
							print('test1')
							tax_name = tax_names[0]
						else:
							print('test2')
							tax_name = process.extract(full_tax_name.replace(' ','_'), tax_names, limit=1)[0]
						print('1')
						print(tax_name)
						tax_name = tax_name[0].replace('[','').replace(']','')
						print('_______________')
						print(tax_name)
						print('_______________')
						if tax_name.replace(' ','_') not in ref_seq.keys():
							print(2)
							print('missing')
							print(3)
							print(tax_name)
							continue
						ref_ac = ref_seq[tax_name.replace(' ','_')][0]
						pos_seq[ref_ac] = sequences[key].replace('-','')

					else:
						new_key = key+';PI='+str(counter/len(ref))
						neg_seq[new_key]=sequences[key]
			# get the ammount of positive and negatives and write to stat file
			pos_count = len(pos_seq)-1
			neg_count = len(neg_seq)
			stat_line = name + '\t'+str(pos_count)+'\t'+str(neg_count)+'\t'+str(count_cor_glut)+'\n'
			stat.write(stat_line)
			# open the file for processing with CopraRNA or IntaRNA
			outfile = open('/data/jannik/data/results/mafft/filtered/'+name,'w')
			for key in pos_seq.keys():
				out_key = '>'+key+'\n'
				outfile.write(out_key)
				outfile.write(pos_seq[key]+'\n')
			outfile.close()
			count_end -=1
			x = [count_end,stat_count-count_end]
			name_file.write(name+'\n')
			name_file.close()
			names.remove(name)
			print(str(count_end)+' of '+str(stat_count)+' to go!')
	except:
		# if an error occures, wait for a time and restart
		time.sleep(60)
		continue

# close the stat file
stat.close()
