# file for downloading raw data from RNAseq : https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6602/E-MTAB-6602/

from ftplib import FTP
import os

ftp = FTP('ftp.sra.ebi.ac.uk')
ftp.login()
ftp.retrlines('LIST')
ftp.cwd('/')

file = open('E-MTAB-6602.sdrf.txt','r')
file.readline()
for line in file:
        ftp_address = line.split('\t')[28].replace('ftp://ftp.sra.ebi.ac.uk/','')
        name = ftp_address[ftp_address.rfind('/')+1:]
        print(name)
        folder = line.split('\t')[26]
        folder = folder[:folder.rfind('_R')]
        if not os.path.exists(folder):
                os.mkdir(folder)
        ftp_address = ftp_address[:ftp_address.rfind('/')]
        ftp.cwd(ftp_address)
        lf = open(folder + '/' + name,'wb')
        ftp.retrbinary('RETR '+name,lf.write, 8*1024)
        lf.close()
        ftp.cwd('/')
        
        
