# script for downloading C. glutamicum ATCC 13032 RefSeq
import Bio
from Bio import Entrez

Entrez.email = 'st145304@stud.uni-stuttgart.de'

handle = Entrez.efetch(db='nuccore',id='62388892',rettype='fasta',retmode='text')

file = open('refseq.fasta', 'w')
file.write(handle.read())
file.close()

