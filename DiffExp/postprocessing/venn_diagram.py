from matplotlib_venn import venn3
from matplotlib_venn import venn2
from matplotlib import pyplot as plt
import numpy as np
import os
# Script for creation of Venn-Diagrams and seperation of the Data from each condition in seperate files
# RNA means the data from Lange et al 2018, sRNA means the data focused on the sRNAs
# NOISeq Results
# Upregulated-RNA
file = open('results/RNA/upregulated.csv','r')

anaerob = []
microaerob1 = []
aerob = []

for entry in file:
    entry = entry.replace('\n','').replace("\"",'').split(',')
    if entry[2] == 'Aerob':
        aerob.append(entry[1])
    elif entry[2] == 'Microaerob':
        microaerob.append(entry[1])
    elif entry[2] == 'Anaerob':
        anaerob.append(entry[1])
file.close()
anaerob = set(anaerob)
microaerob = set(microaerob)
aerob = set(aerob)

all = set()
for entry in aerob:
    all.add(entry)
for entry in anaerob:
    all.add(entry)
for entry in microaerob:
    all.add(entry)

file = open('results/RNA/aerob_up.txt','w')
for entry in aerob:
    file.write(entry+'\n')
file.close()

file = open('results/RNA/anaerob_up.txt','w')
for entry in anaerob:
    file.write(entry+'\n')
file.close()

file = open('results/RNA/microaerob_up.txt','w')
for entry in microaerob:
    file.write(entry+'\n')
file.close()


#
#   sRNA-Upregulated
#
file = open('results/sRNA/upregulated.csv','r')

sanaerob = []
smicroaerob = []
saerob = []
for entry in file:
    entry = entry.replace('\n','').replace("\"",'').split(',')
    if entry[2] == 'Aerob':
        saerob.append(entry[1])
    elif entry[2] == 'Microaerob':
        smicroaerob.append(entry[1])
    elif entry[2] == 'Anaerob':
        sanaerob.append(entry[1])
file.close()
sanaerob = set(sanaerob)
smicroaerob = set(smicroaerob)
saerob = set(saerob)



sall = set()
for entry in saerob:
    sall.add(entry)
for entry in sanaerob:
    sall.add(entry)
for entry in smicroaerob:
    sall.add(entry)

file = open('results/sRNA/aerob_up.txt','w')
for entry in saerob:
    file.write(entry+'\n')
file.close()

file = open('results/sRNA/anaerob_up.txt','w')
for entry in sanaerob:
    file.write(entry+'\n')
file.close()

file = open('results/sRNA/microaerob_up.txt','w')
for entry in smicroaerob:
    file.write(entry+'\n')
file.close()

fig = plt.figure(figsize=(5,10), dpi = 300)
plt.rcParams.update({'font.size': 20})
ax1 = fig.add_subplot(211, title = 'Lange et al \'18')
plt.rcParams.update({'font.size': 14})
venn = venn3([aerob,anaerob,microaerob], ('Aerob', 'Anaerob', 'Microaerob'),ax = ax1)
plt.rcParams.update({'font.size': 20})
ax2 = fig.add_subplot(212, title = 'small RNA Focus')
plt.rcParams.update({'font.size': 14})
venn = venn3([saerob,sanaerob,smicroaerob], ('Aerob', 'Anaerob', 'Microaerob'),ax =ax2)
fig.savefig('results/Upregulated.png')
plt.close(fig)


# Venndiagramm of all in RNA and sRNA focus - upregulated

fig = plt.figure(figsize=(5,5), dpi = 300)
ax1 = fig.add_subplot(111, title = 'Differentially expressed RNAs')
venn = venn2([all,sall], ('RNA','sRNA'),ax =ax1)
fig.savefig('results/up.png')
plt.close(fig)


# Downregulated
file = open('results/RNA/downregulated.csv','r')

anaerob = []
microaerob = []
aerob = []

for entry in file:
    entry = entry.replace('\n','').replace("\"",'').split(',')
    if entry[2] == 'Aerob':
        aerob.append(entry[1])
    elif entry[2] == 'Microaerob':
        microaerob.append(entry[1])
    elif entry[2] == 'Anaerob':
        anaerob.append(entry[1])
file.close()

anaerob = set(anaerob)
microaerob = set(microaerob)
aerob = set(aerob)

all = set()
for entry in aerob:
    all.add(entry)
for entry in anaerob:
    all.add(entry)
for entry in microaerob:
    all.add(entry)
file.close()

file = open('results/RNA/aerob_down.txt','w')
for entry in aerob:
    file.write(entry+'\n')
file.close()

file = open('results/RNA/anaerob_down.txt','w')
for entry in anaerob:
    file.write(entry+'\n')
file.close()

file = open('results/RNA/microaerob_down.txt','w')
for entry in microaerob:
    file.write(entry+'\n')
file.close()




#
#   sRNA-Downregulated
#
file = open('results/sRNA/downregulated.csv','r')

sanaerob = []
smicroaerob = []
saerob = []
for entry in file:
    entry = entry.replace('\n','').replace("\"",'').split(',')
    if entry[2] == 'Aerob':
        saerob.append(entry[1])
    elif entry[2] == 'Microaerob':
        smicroaerob.append(entry[1])
    elif entry[2] == 'Anaerob':
        sanaerob.append(entry[1])
file.close()
len(sanaerob)
sanaerob = set(sanaerob)
smicroaerob = set(smicroaerob)
saerob = set(saerob)

sall = set()
for entry in saerob:
    sall.add(entry)
for entry in sanaerob:
    sall.add(entry)
for entry in smicroaerob:
    sall.add(entry)

file = open('results/sRNA/aerob_down.txt','w')
for entry in saerob:
    file.write(entry+'\n')
file.close()

file = open('results/sRNA/anaerob_down.txt','w')
for entry in sanaerob:
    file.write(entry+'\n')
file.close()

file = open('results/sRNA/microaerob_down.txt','w')
for entry in smicroaerob:
    file.write(entry+'\n')
file.close()

fig = plt.figure(figsize=(5,10), dpi = 300)
plt.rcParams.update({'font.size': 20})
ax1 = fig.add_subplot(211, title = 'Lange et al \'18')
plt.rcParams.update({'font.size': 14})
venn = venn3([aerob,anaerob,microaerob], ('Aerob', 'Anaerob', 'Microaerob'),ax =ax1)
plt.rcParams.update({'font.size': 20})
ax2 = fig.add_subplot(212, title = 'small RNA Focus')
plt.rcParams.update({'font.size': 14})
venn = venn3([saerob,sanaerob,smicroaerob], ('Aerob', 'Anaerob', 'Microaerob'),ax =ax2)
fig.savefig('results/Downregulated.png')
plt.close(fig)

# Venndiagramm of all in RNA and sRNA focus - upregulated

fig = plt.figure(figsize=(5,5), dpi = 300)
ax1 = fig.add_subplot(111, title = 'Differentially expressed RNAs')
venn = venn2([all,sall], ('RNA','sRNA'),ax =ax1)
fig.savefig('results/down.png')
plt.close(fig)

# Overall Differentially expressed
file = open('results/RNA/all.csv','r')

anaerob = []
microaerob = []
aerob = []

for entry in file:
    entry = entry.replace('\n','').replace("\"",'').split(',')
    if entry[2] == 'Aerob':
        aerob.append(entry[1])
    elif entry[2] == 'Microaerob':
        microaerob.append(entry[1])
    elif entry[2] == 'Anaerob':
        anaerob.append(entry[1])
file.close()
anaerob = set(anaerob)
microaerob = set(microaerob)
aerob = set(aerob)


all = set()
for entry in aerob:
    all.add(entry)
for entry in anaerob:
    all.add(entry)
for entry in microaerob:
    all.add(entry)

# Overall Differentially expressed-sRNA
file = open('results/sRNA/all.csv','r')

sanaerob = []
smicroaerob = []
saerob = []
for entry in file:
    entry = entry.replace('\n','').replace("\"",'').split(',')
    if entry[2] == 'Aerob':
        saerob.append(entry[1])
    elif entry[2] == 'Microaerob':
        smicroaerob.append(entry[1])
    elif entry[2] == 'Anaerob':
        sanaerob.append(entry[1])
file.close()
len(sanaerob)
sanaerob = set(sanaerob)
smicroaerob = set(smicroaerob)
saerob = set(saerob)



sall = set()
for entry in saerob:
    sall.add(entry)
for entry in sanaerob:
    sall.add(entry)
for entry in smicroaerob:
    sall.add(entry)
len(sall)

fig = plt.figure(figsize=(5,10), dpi = 300)
plt.rcParams.update({'font.size': 20})
ax1 = fig.add_subplot(211, title = 'Lange et al \'18')
plt.rcParams.update({'font.size': 14})
venn = venn3([aerob,anaerob,microaerob], ('Aerob', 'Anaerob', 'Microaerob'),ax =ax1)
plt.rcParams.update({'font.size': 20})
ax2 = fig.add_subplot(212, title = 'small RNA Focus')
plt.rcParams.update({'font.size': 14})
venn = venn3([saerob,sanaerob,smicroaerob], ('Aerob', 'Anaerob', 'Microaerob'),ax =ax2)
fig.savefig('results/all.png')
plt.close(fig)

fig = plt.figure(figsize=(5,5), dpi = 300)
ax1 = fig.add_subplot(111, title = 'Differentially expressed RNAs')
venn = venn2([all,sall], ('RNA','sRNA'),ax =ax1)
fig.savefig('results/all_two_cond.png')
plt.close(fig)
