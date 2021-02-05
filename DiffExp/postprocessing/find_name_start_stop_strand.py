    # script for extracting name start stop strand and featuretype from featureCount
    # Countdata

    def find_asRNA(name):
        file1 = open('counts/RNA/counts_asRNAs.txt','r')
        file2 = open('counts/sRNA/counts_asRNAs.txt', 'r')

        genid1 = []
        genid2 = []

        start1 = []
        start2 = []

        stop1 = []
        stop2 = []

        strand1 = []
        strand2 = []

        for entry in file1:
            entry = entry.replace('\n','').split('\t')
            if entry[0][0] != '#':

                genid1.append(entry[0])
                start1.append(entry[2])
                stop1.append(entry[3])
                strand1.append(entry[4])

        for entry in file2:
            entry = entry.replace('\n','').split('\t')
            if entry[0][0] != '#':

                genid1.append(entry[0])
                start1.append(entry[2])
                stop1.append(entry[3])
                strand1.append(entry[4])

        if name in genid1:
            position1 = genid1.index(name)
        if name in genid2:
            position2 = genid2.index(name)

        if 'position1' in locals():
            # list to return
            ret_list = [genid1[position1],start1[position1],stop1[position1], strand1[position1], 'asRNA']
            return ret_list
        elif 'position2' in locals():
            ret_list = [genid2[position2],start2[position2],stop2[position2], strand2[position2], 'asRNA']
            return ret_list
        else:
            return 'NA'


    def find_mRNA(name):
        file1 = open('counts/RNA/counts_mRNAs.txt','r')
        file2 = open('counts/sRNA/counts_mRNAs.txt', 'r')

        genid1 = []
        genid2 = []

        start1 = []
        start2 = []

        stop1 = []
        stop2 = []

        strand1 = []
        strand2 = []

        for entry in file1:
            entry = entry.replace('\n','').split('\t')
            if entry[0][0] != '#':

                genid1.append(entry[0])
                start1.append(entry[2])
                stop1.append(entry[3])
                strand1.append(entry[4])

        for entry in file2:
            entry = entry.replace('\n','').split('\t')
            if entry[0][0] != '#':

                genid1.append(entry[0])
                start1.append(entry[2])
                stop1.append(entry[3])
                strand1.append(entry[4])
        if name in genid1:
            position1 = genid1.index(name)
        if name in genid2:
            position2 = genid2.index(name)

        if 'position1' in locals():
            # list to return
            ret_list = [genid1[position1],start1[position1],stop1[position1], strand1[position1], 'mRNA']
            return ret_list
        elif 'position2' in locals():
            ret_list = [genid2[position2],start2[position2],stop2[position2], strand2[position2], 'mRNA']
            return ret_list
        else:
            return 'NA'


    def find_ncRNA(name):
        file1 = open('counts/RNA/counts_ncRNAs.txt','r')
        file2 = open('counts/sRNA/counts_ncRNAs.txt', 'r')

        genid1 = []
        genid2 = []

        start1 = []
        start2 = []

        stop1 = []
        stop2 = []

        strand1 = []
        strand2 = []

        for entry in file1:
            entry = entry.replace('\n','').split('\t')
            if entry[0][0] != '#':

                genid1.append(entry[0])
                start1.append(entry[2])
                stop1.append(entry[3])
                strand1.append(entry[4])

        for entry in file2:
            entry = entry.replace('\n','').split('\t')
            if entry[0][0] != '#':

                genid1.append(entry[0])
                start1.append(entry[2])
                stop1.append(entry[3])
                strand1.append(entry[4])

        if name in genid1:
            position1 = genid1.index(name)
        if name in genid2:
            position2 = genid2.index(name)

        if 'position1' in locals():
            # list to return
            ret_list = [genid1[position1],start1[position1],stop1[position1], strand1[position1], 'ncRNA']
            return ret_list
        elif 'position2' in locals():
            ret_list = [genid2[position2],start2[position2],stop2[position2], strand2[position2], 'ncRNA']
            return ret_list
        else:
            return 'NA'


    def extract_features(path):
        file = open(path, 'r')

        features = []
        for entry in file:
            entry = entry.replace('\n','').split(',')
            ret_list1 = find_asRNA(entry[1].replace('\"',''))
            ret_list2 = find_mRNA(entry[1].replace('\"',''))
            ret_list3 = find_ncRNA(entry[1].replace('\"',''))
            if ret_list1 == 'NA':
                ret_list1 = []
            if ret_list2 == 'NA':
                ret_list2 = []
            if ret_list3 == 'NA':
                ret_list3 = []
            if ret_list1 and not ret_list2 and not ret_list3:
                ret_list1.append(entry[2].replace('\"',''))
                features.append(ret_list1)
            elif ret_list1 and ret_list2 and not ret_list3:
                print('asRNA','mRNA')
            elif ret_list1 and ret_list2 and ret_list3:
                print('asRNA','mRNA', 'ncRNA')
            elif ret_list1 and not ret_list2 and ret_list3:
                print('asRNA', 'ncRNA')
            elif not ret_list1 and ret_list2 and not ret_list3:
                ret_list2.append(entry[2].replace('\"',''))
                features.append(ret_list2)
            elif not ret_list1 and ret_list2 and ret_list3:
                print('mRNA','ncRNA')
            elif not ret_list1 and not ret_list2 and ret_list3:
                ret_list3.append(entry[2].replace('\"',''))
                features.append(ret_list3)
        file.close()
        outpath = path.replace('.csv','_feature.txt')
        file = open(outpath,'w')
        print('EXITING')
        print(features)
        for line in features:
            file.write(line[0]+'\t'+line[1]+'\t'+line[2]+'\t'+line[3]+'\t'+line[4]+'\t'+line[5]+'\n')
        file.close()

    path_list = ['results/sRNA/downregulated_tp.csv','results/RNA/downregulated_tp.csv',
     'results/sRNA/upregulated_tp.csv', 'results/RNA/upregulated_tp.csv']

    for path in path_list:
        extract_features(path)
