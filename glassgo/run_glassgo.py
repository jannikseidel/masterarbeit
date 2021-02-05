# Script for running glassgo with the noncoding RNAs
import pandas as pd
from Bio.Seq import Seq
import subprocess as sp
import os
import time
# Set working directory to /
os.chdir('/')
# Specify the path to the lists with the RNAs in Format id \t start \t stop \t strand \t featuretype \n
path_list = ['/data/jannik/data/results/diffexp/all_merged_feature_nc.txt']


# prepare list of lists for generation of pandas dataframe
d_frame_prepare = []
for path in path_list:
    file = open(path,'r')
    for entry in file:
        d_frame_prepare.append(entry.replace('\n','').split('\t'))
    file.close()
# Create dataframe and drop duplicates
d_frame = pd.DataFrame(d_frame_prepare,columns= ['id','start','stop','strand','featuretype','timepoint','direction','log2fold change'])
d_frame = d_frame.drop('timepoint', axis=1)
d_frame = d_frame.drop_duplicates(subset = 'id')
print(d_frame)
# path to reference genome
ref_path = '/data/jannik/data/sequences/refseq/refseq.fasta'

file = open(ref_path,'r')
ref_seq = []
seq = ''
for entry in file:
    if entry[0] == '>':
        ref_seq.append(entry.replace('\n',''))
    elif entry[0] != '>' and entry[0] != '':
        seq += entry.replace('\n','')
file.close()
seq = Seq(seq)
ref_seq.append(seq)
counter = 0
start_time = time.perf_counter()

list_names = list(d_frame.index.values)
# run glassgo for the ncRNAs
while list_names:
    try:
        for entry in list_names:
            feature_object = d_frame.loc[entry,:]
            print(feature_object)
            command = 'docker exec -u root -it 0cc3f9db46300625278a77c98c9677222b7f5ddd8292d218af0998694a77c95e ./GLASSgo.py -t 26 -d /BLAST_NT/nt -i /SRNA_INPUT/sequence.fasta -o /results/'
            # writing the DNA Sequence of the nc RNA to fasta file to be processed by
            # glassgo and run GLASSgo

            if feature_object['strand'] == '+':
                seq = ref_seq[1][int(feature_object['start'])-1:int(feature_object['stop'])]
                file = open('/data/jannik/data/results/diffexp/sequence.fasta','w')
                file.write(ref_seq[0]+';'+feature_object['id']+'\n')
                file.write(str(seq)+'\n')
                file.close()
                command = command + feature_object['id']+'.fasta'
                command = command.split(' ')
                sp.call(command)
            elif feature_object['strand'] == '-':
                seq = ref_seq[1][int(feature_object['start'])-1:int(feature_object['stop'])].reverse_complement()
                file = open('/data/jannik/data/results/diffexp/sequence.fasta','w')
                file.write(ref_seq[0]+';'+feature_object['id']+'\n')
                file.write(str(seq)+'\n')
                file.close()
                command = command + feature_object['id']+'.fasta'
                command = command.split(' ')
                sp.call(command)

            counter += 1
            print('______________________________')
            print('GLASSgo - time taken: ' + str((time.perf_counter()-start_time)/3600) + ' h est. to go:' + str((time.perf_counter()-start_time)/(counter*3600)*(len(d_frame.index.values)-counter))+ 'h')
            print('----------',counter,'---------')
            list_names.remove(entry)
    except:
        continue
print('finished GLASSgo')
