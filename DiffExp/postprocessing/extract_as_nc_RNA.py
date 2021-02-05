# script for extraction of nc and asRNA from file with annotation

path_list = ['results/sRNA/downregulated_tp_feature.txt','results/RNA/downregulated_tp_feature.txt',
 'results/sRNA/upregulated_tp_feature.txt', 'results/RNA/upregulated_tp_feature.txt']

def extract(path):
    file = open(path)
    as_nc = []
    for line in file:
        line = line.replace('\n','').split('\t')
        if line[4] == 'asRNA' or line[4] == 'ncRNA':
            as_nc.append(line)
    file.close()
    out_path = path.replace('.txt','_as_nc.txt')
    file = open(out_path,'w')
    for line in as_nc:
        file.write(line[0]+'\t'+line[1]+'\t'+line[2]+'\t'+line[3]+'\t'+line[4]+'\n')
    file.close()

for entry in path_list:
    extract(entry)

def extract_nc(path):
    file = open(path)
    as_nc = []
    for line in file:
        line = line.replace('\n','').split('\t')
        if line[4] == 'ncRNA':
            as_nc.append(line)
    file.close()
    out_path = path.replace('.txt','_nc.txt')
    file = open(out_path,'w')
    for line in as_nc:
        file.write(line[0]+'\t'+line[1]+'\t'+line[2]+'\t'+line[3]+'\t'+line[4]+'\t'+line[5]+'\n')
    file.close()

for entry in path_list:
    extract_nc(entry)
