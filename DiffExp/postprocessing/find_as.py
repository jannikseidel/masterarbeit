# script for finding the gene corresponding to an asRNA
import Bio
import pandas as pd
import math
# extracting the features out of the gff file
annot = open('NC_006958.1.gff3','r')
name = []
start = []
stop = []
strand = []
species = []
locus_tag = []

for line in annot:
    line = line.strip('\n').split('\t')
    if line[0][0] == '#':
        continue

    if line[1] == 'RefSeq':
        id = line[8][line[8].find('=')+1:line[8].find(';')]
        name.append(id)
        start.append(int(line[3]))
        stop.append(int(line[4]))
        strand.append(line[6])
        species.append(line[2])
        if 'locus_tag' in line[8]:
            tag_start = line[8].find('locus_tag=')+len('locus_tag=')
            tag_stop = line[8].rfind(';',tag_start)
            print(tag_start)
            print(line[8])
            print(tag_stop)
            tag = line[8][tag_start:tag_stop]
            print(tag)
            locus_tag.append(tag)
        else:
            locus_tag.append('NA')
annot.close()
# create dataframe for gff annotation file
gff = pd.DataFrame({'genid' : name, 'start':start,'stop':stop, 'strand':strand, 'species':species,'locus_tag':locus_tag})
gff = gff.sort_values(by='start')

def find_start(value,annot_start):
    # find position which is potential start of antisense RNA
    searching = True
    position = 0
    position1 = len(annot_start)-1
    check = 0
    while searching :
        check = position+position1
        start_search = math.floor((position1+position)/2)
        if int(annot_start.values[start_search]) > value and start_search > 0:
            position1 = start_search -1
        elif int(annot_start.values[start_search]) < value and start_search < len(annot_start):
            position = start_search +1

        if check == (position + position1):

            if int(annot_start.values[position]) <= value:
                return position
            elif int(annot_start.values[position1]) <= value:
                return position1

def dest_as_trans(start,stop, annot_start,annot_stop):
    # destinguish a small RNA between cis and trans -> What if same strand! [XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX]
    pos_start = find_start(start,annot_start)

    if annot_stop.values[pos_start] >= stop:
        anti = 'anti'

        return pos_start,anti
    elif annot_start.values[pos_start] <= start and annot_stop.values[pos_start] > start and annot_stop.values[pos_start] <= stop:
        anti = 'anti'

        return pos_start,anti
    elif annot_start.values[pos_start] >= start and annot_start.values[pos_start] < stop:
        anti = 'anti'

        return pos_start,anti
    elif annot_stop.values[pos_start] < start:
        trans = 'trans'

        return '',trans

path_list = ["results/sRNA/downregulated_tp_feature.txt","results/RNA/downregulated_tp_feature.txt", "results/sRNA/upregulated_tp_feature.txt","results/RNA/upregulated_tp_feature.txt"]

for path in path_list:
    antis = []
    trans = []
    file = open(path,'r')
    for line in file:
        line = line.replace('\n','').split('\t')
        if line[4] == 'asRNA':

            line1 = int(line[1])
            line2 = int(line[2])
            start = gff['start']
            stop = gff['stop']
            pos,state = dest_as_trans(line1,line2,start,stop)
            if pos == '':
                trans.append(line)
            else:
                if line[0] != gff['genid'].values[pos]:
                    line[4] = line[4]+';antisenseto='+gff['genid'].values[pos]+';locus_tag='+gff['locus_tag'].values[pos]

                    antis.append(line)
    file.close()
    out_path = path.replace('.txt','_as.txt')
    file = open(out_path,'w')
    for line in antis:

        file.write(line[0]+'\t'+line[1]+'\t'+line[2]+'\t'+line[3]+'\t'+line[4]+'\t'+line[5]+'\n')
    file.close()
