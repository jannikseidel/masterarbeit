import Bio.Seq as Seq

# open reference genome
file = open('/home/janniksl/Studium/master/Masterarbeit/sequences/refseq/refseq.fasta')

seq = ['','']
for entry in file:
    if entry[0] == '>':
        seq[0] = entry.replace('\n','')
    else:
        seq[1] += entry.replace('\n','')
seq[1] = Seq.Seq(seq[1])

print(seq[1][3135511:3135546].reverse_complement())
