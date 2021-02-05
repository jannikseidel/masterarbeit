# script for examining the filtered data and moving them to the corresponding
# subfolders according to the ammount of sequences with a PI >= 80% compared
# to the sequence of the organism of interest

import pandas as pd
import os
import shutil

# opening the csv with the stats on the MAFFT MSAs

file = pd.read_csv('/data/jannik/data/results/mafft/filtered/stats.txt',sep='\t',lineterminator='\n')

counter_pos = 0
counter_neg = 0
os.chdir('/')

# itterating over the entries in the stats file and trying to move the files to
# a corresponding subfolder (classified upon the number of sequences >= 80% PI,
# # >= 3 -> CopraRNA, # < 3 IntaRNA)

for i in range(0,len(file)):
    row = file.iloc[[i]]
    if int(row['positives']) >=3:
        try:
            os.chdir('/data/jannik/data/results/mafft/filtered/CopraRNA/')
            counter_pos += 1
            destination = '/data/jannik/data/results/mafft/filtered/CopraRNA/'
            source = '/data/jannik/data/results/mafft/filtered/'
            name = str(row.iloc[0]['#filename'])

            try:
                os.mkdir(name.replace('.fasta',''), mode=0o777)
            except:
                print('Examining the filtered Sequences-folder not created')
            destination = destination + name.replace('.fasta','') + '/' + name
            source += name
            shutil.copy(source,destination)
        except:
            continue
    else:
        try:
            os.chdir('/data/jannik/data/results/mafft/filtered/IntaRNA/')

            destination = '/data/jannik/data/results/mafft/filtered/IntaRNA/'
            source= '/data/jannik/data/results/mafft/filtered/'
            name = str(row.iloc[0]['#filename'])

            try:
                os.mkdir(name.replace('.fasta',''), mode=0o777)
            except:
                print('Examining the filtered Sequences-folder not created')
            destination = destination +name.replace('.fasta','') + '/' + name
            source += name
            shutil.copy(source,destination)
            counter_neg += 1
        except:
            continue
