file = open('mrna.fasta','r')
seq = ''
for entry in file:
    if entry[0] == '>':
        continue
    else:
        seq += entry.replace('\n','')
print(len(seq))
