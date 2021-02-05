# script for extracting the genes annotated in a gff from ref_seq genome to
# a file only containing the mRNA-Sequences to be processed by IntaRNA
import Bio.Seq as Seq


# loading the relevant data from gff into list with [[start,stop,strand]]
file = open('/home/janniksl/Studium/master/Masterarbeit/annotations/NC_006958.1.gff3', 'r')

# counter for counting the lines of file already processed
counter = 0

# extracting the relevant annotaions to a list
annotation = []

for entry in file:
	entry = entry.replace('\n', '').split('\t')
	if counter >= 4 and len(entry) > 1:
		if entry[1] == 'RefSeq':
			name = entry[8].replace('ID=','')
			name = name.split(';')[1]
			name = name.replace('Name=','')
			#print(entry)
			annotation.append([entry[3],entry[4],entry[6], name])
	counter += 1

# open reference genome
file = open('/home/janniksl/Studium/master/Masterarbeit/sequences/refseq/refseq.fasta')

seq = ['','']
for entry in file:
    if entry[0] == '>':
        seq[0] = entry.replace('\n','')
    else:
        seq[1] += entry.replace('\n','')
seq[1] = Seq.Seq(seq[1])

# transcribe the annotated features from dna to mrna
mrna = []
for entry in annotation:
	if entry[2] == '+':
		start = int(entry[0])
		stop = int(entry[1])
		mr = seq[1][start:stop].complement().transcribe()
		mrna.append([mr,entry[3]])
	if entry[2] == '-':
		start = int(entry[0])
		stop = int(entry[1])
		mr = seq[1][start:stop][::-1].transcribe()
		mrna.append([mr,entry[3]])

# open the files where the new sequence is stored and the corresponding annotation
file = open('mrna.fasta', 'w')
new_file = open('mrna_annot.txt','w')
file.write('>C. glutamicum-transcribed \n')
len_0 = 1
len_1 = 0
seq = ''
for entry in mrna:

	len_1 = len(entry[0])+ len_0-1
	pos = str(len_0)+'\t'+str(len_1)
	len_0 = len_1+1
	seq += str(entry[0])
	print(pos)
	print(entry[1])
	print(len_1,len(seq))
	new_file.write(pos+'\t'+entry[1]+'\n')
	file.write(str(entry[0]))

file.close()

file = open('mrna_annot.txt','r')

for entry in file:
	print(entry)
428488	
