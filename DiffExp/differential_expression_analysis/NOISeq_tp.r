##
library('NOISeq')
library("gplots")
library(reticulate)

#####################################################################################################
# sRNA - Grouped
### SETTINGS
# set (full path) the directory containing the counts (without dash at the end)
dir = "C:/Users/seide/Dropbox/Studium/Master/Masterarbeit/c_glutamicum_atcc_13032/DiffExp/counts/sRNA" 
nr_datasets = 6 # number of datasets


results_dir = "C:/Users/seide/Dropbox/Studium/Master/Masterarbeit/c_glutamicum_atcc_13032/DiffExp/results/sRNA"
setwd(results_dir)


### DATA PREPARATION
fnames <- list.files(path=dir,full.names=T,include.dirs = FALSE) # retrieve abs path of count files
fnames <- fnames[1:3]
data <- lapply(fnames,read.table, sep="\t", header=TRUE) # read count files

types <- basename(fnames)
types <- c("asRNAs","mRNAs","ncRNAs")
types <- rep(types,sapply(data,nrow))


# merge count files (list of lists) into one
data <- do.call("rbind",data)

mylength <- data.frame(data$Geneid,data$Length)

counts <- data.frame(data$X.data.jannik.data.rnaseq.sRNA_focus.JL1_sRNA_3h_AGTCAA_L001_R1_001.mapped.sorted.bam,data$X.data.jannik.data.rnaseq.sRNA_focus.JL2_sRNA_5h_AGTTCC_L001_R1_001.mapped.sorted.bam,data$X.data.jannik.data.rnaseq.sRNA_focus.JL3_sRNA_7h_ATGTCA_L001_R1_001.mapped.sorted.bam,data$X.data.jannik.data.rnaseq.sRNA_focus.JL4_sRNA_9h_CCGTCC_L001_R1_001.mapped.sorted.bam,data$X.data.jannik.data.rnaseq.sRNA_focus.JL5_sRNA_11h_GTAGAG_L001_R1_001.mapped.sorted.bam,data$X.data.jannik.data.rnaseq.sRNA_focus.JL6_sRNA_13h_GTCCGC_L001_R1_001.mapped.sorted.bam)
rownames(counts) <-data$Geneid
col_names <- c('3h','5h','7h','9h','11h','13h')
colnames(counts) <- col_names

mybiotypes <- data.frame(types)
rownames(mybiotypes) <- data$Geneid

myfactors <- data.frame(Organism = c('c_glutamicum','c_glutamicum','c_glutamicum','c_glutamicum','c_glutamicum','c_glutamicum'),Group = c('aerob-Reference', 'aerob', 'microaerob', 'microaerob', 'microaerob', 'anaerob'))

mydata <- readData(data = counts, length = mylength, factors = myfactors, biotype = mybiotypes)


# Quality Controll
mycountsbio = dat(mydata, factor = NULL, type = "countsbio")
mycountsbio = dat(mydata, factor = NULL, type = "countsbio")
mysaturation <- dat(mydata, k = 0, ndepth = 7, type = "saturation")
mylengthbias <- dat(mydata, factor = "Group", type = "lengthbias")
myPCA <- dat(mydata, type = "PCA")

png('explorative.png',width=10, height = 27, units = 'cm',  res= 300)
nf <- layout(matrix(c(1,2,3,4),ncol=1), widths=c(4,4,4,4), heights=c(4,4,4,4), TRUE)
explo.plot(mysaturation, toplot = 1, samples = 1:6, yleftlim = NULL, yrightlim = NULL)
explo.plot(mylengthbias, samples = NULL, toplot = "global")
explo.plot(myPCA, factor = "Group")
explo.plot(mycountsbio, toplot = 1, samples = NULL, plottype = "boxplot")
dev.off()

#Steps included in noiseq function
# Normalisation with TMM

#myTMM <- tmm(assayData(mydata)$exprs,long = mylength, lc = 1)

# filtering of Genes with low counts -> noise reduction
#myfilt <- filtered.data(counts, factor = myfactors$Group, norm = FALSE, depth = NULL, method = 1, cv.cutoff = 100, cpm = 1,p.adj = 'fdr')
prob_threshold <- 0.92
# Calculating results aerob vs aerob
myresults_a_a <- noiseq(mydata, factor = 'Group', conditions = c('aerob-Reference', 'aerob'), k = NULL, norm = 'tmm',pnr = 0.2, v = 0.02, lc = 1, replicates = 'no')
myresults_a_a.deg <- degenes(myresults_a_a, q = prob_threshold, M = NULL)
myresults_a_a.deg1 <- degenes(myresults_a_a, q = prob_threshold, M = 'up')
myresults_a_a.deg2 <- degenes(myresults_a_a, q = prob_threshold, M = 'down')


# Calculating results aerob vs anaerob
myresults_a_an <- noiseq(mydata, factor = 'Group', conditions = c('aerob-Reference', 'anaerob'), k = NULL, norm = 'tmm',pnr = 0.2, v = 0.02, lc = 1, replicates = 'no')
myresults_a_an.deg <- degenes(myresults_a_an, q = prob_threshold, M = NULL)
myresults_a_an.deg1 <- degenes(myresults_a_an, q = prob_threshold, M = 'up')
myresults_a_an.deg2 <- degenes(myresults_a_an, q = prob_threshold, M = 'down')


# Calculating results aerob vs microaerob
myresults_a_mic <- noiseq(mydata, factor = 'Group', conditions = c('aerob-Reference', 'microaerob'), k = NULL, norm = 'tmm',pnr = 0.2, v = 0.02, lc = 1, replicates = 'no')
myresults_a_mic.deg <- degenes(myresults_a_mic, q = prob_threshold, M = NULL)
myresults_a_mic.deg1 <- degenes(myresults_a_mic, q = prob_threshold, M = 'up')
myresults_a_mic.deg2 <- degenes(myresults_a_mic, q = prob_threshold, M = 'down')


# Calculating results anaerob vs microaerob
#myresults_an_mic <- noiseq(mydata, factor = 'Group', conditions = c('microaerob', 'anaerob'), k = NULL, norm = 'tmm',pnr = 0.2, v = 0.02, lc = 1, replicates = 'no')
#myresults_an_mic.deg <- degenes(myresults_an_mic, q = prob_threshold, M = NULL)
#myresults_an_mic.deg1 <- degenes(myresults_an_mic, q = prob_threshold, M = 'up')
#myresults_an_mic.deg2 <- degenes(myresults_an_mic, q = prob_threshold, M = 'down')

# Plotting of the MD DIstribution with highlighting of the DE genes
png('MD.png',width=10, height = 27, units = 'cm',  res= 300)
nf <- layout(matrix(c(1,2,3,4),ncol=1), widths=c(4,4,4,4), heights=c(4,4,4,4), TRUE)
DE.plot(myresults_a_a, q = prob_threshold, graphic = 'MD')
DE.plot(myresults_a_mic, q = prob_threshold, graphic = 'MD')
DE.plot(myresults_a_an, q = prob_threshold, graphic = 'MD')
#DE.plot(myresults_an_mic, q = prob_threshold, graphic = 'MD')
dev.off()
png('Expression.png',width=15, height = 27, units = 'cm',  res= 300)
nf<-layout(matrix(c(1,2,3,4),ncol=1), widths=c(4,4,4,4), heights=c(4,4,4,4), TRUE)
DE.plot(myresults_a_a, q = prob_threshold, graphic = "expr", log.scale = TRUE)
DE.plot(myresults_a_an, q = prob_threshold, graphic = "expr", log.scale = TRUE)
DE.plot(myresults_a_mic, q = prob_threshold, graphic = "expr", log.scale = TRUE)
#DE.plot(myresults_an_mic, q = prob_threshold, graphic = "expr", log.scale = TRUE)
dev.off()
png('Distribution.png',width=10, height = 27, units = 'cm',  res= 300)
nf<-layout(matrix(c(1,2,3,4),ncol=1), widths=c(4,4,4,4), heights=c(4,4,4,4), TRUE)
DE.plot(myresults_a_a, chromosomes = NULL, q = prob_threshold, graphic = "distr")
DE.plot(myresults_a_an, chromosomes = NULL, q = prob_threshold, graphic = "distr")
DE.plot(myresults_a_mic, chromosomes = NULL, q = prob_threshold, graphic = "distr")
#DE.plot(myresults_an_mic, chromosomes = NULL, q = prob_threshold, graphic = "distr")
dev.off()
# extract the feature names of the differential expressed features between the conditions
dif_exp_a_a <- rownames(myresults_a_a.deg)
dif_exp_a_an <- rownames(myresults_a_an.deg)
dif_exp_a_mic <- rownames(myresults_a_mic.deg)
#dif_exp_an_mic <- rownames(myresults_an_mic.deg)

# extract the feature names of the differential expressed features between the conditions
dif_exp_a_a_up <- rownames(myresults_a_a.deg1)
dif_exp_a_an_up <- rownames(myresults_a_an.deg1)
dif_exp_a_mic_up <- rownames(myresults_a_mic.deg1)
#dif_exp_an_mic_up <- rownames(myresults_an_mic.deg1)

# extract the feature names of the differential expressed features between the conditions
dif_exp_a_a_down <- rownames(myresults_a_a.deg2)
dif_exp_a_an_down <- rownames(myresults_a_an.deg2)
dif_exp_a_mic_down <- rownames(myresults_a_mic.deg2)
#dif_exp_an_mic_down <- rownames(myresults_an_mic.deg2)



dif_exp_all <- unique(c(dif_exp_a_a,dif_exp_a_an,dif_exp_a_mic))
dif_exp_all_up <- unique(c(dif_exp_a_a_up,dif_exp_a_an_up,dif_exp_a_mic_up))
dif_exp_all_down <- unique(c(dif_exp_a_a_down,dif_exp_a_an_down,dif_exp_a_mic_down))



length_df <- as.numeric(max(length(dif_exp_a_a_down),length(c(dif_exp_a_mic_down)),length(c(dif_exp_a_an_down))))
a_down <- dif_exp_a_a_down
mic_down <- c(dif_exp_a_mic_down)
an_down <- c(dif_exp_a_an_down)


d_frame <- matrix(nrow = (length(a_down)+length(an_down)+length(mic_down)),ncol = 2) 
counter <- 1
for(entry in a_down){
  d_frame[counter,1] <- entry
  d_frame[counter,2] <- 'Aerob'
  counter <- counter +1
}
for(entry in an_down){
  d_frame[counter,1] <- entry
  d_frame[counter,2] <- 'Anaerob'
  counter <- counter +1
}
for(entry in mic_down){
  d_frame[counter,1] <- entry
  d_frame[counter,2] <- 'Microaerob'
  counter <- counter +1
}


write.csv(d_frame,file = 'downregulated.csv')
  


# Upregulated 
a_up <- dif_exp_a_a_up
mic_up <- dif_exp_a_mic_up
an_up <- dif_exp_a_an_up

d_frame <- matrix(nrow = (length(a_up)+length(an_up)+length(mic_up)),ncol = 2) 
counter <- 1
for(entry in a_up){
  d_frame[counter,1] <- entry
  d_frame[counter,2] <- 'Aerob'
  counter <- counter +1
}
for(entry in an_up){
  d_frame[counter,1] <- entry
  d_frame[counter,2] <- 'Anaerob'
  counter <- counter +1
}
for(entry in mic_up){
  d_frame[counter,1] <- entry
  d_frame[counter,2] <- 'Microaerob'
  counter <- counter +1
}


write.csv(d_frame,file = 'upregulated.csv')

# Upregulated & Downregulated
a_all <- dif_exp_a_a
mic_all <- dif_exp_a_mic
an_all <- dif_exp_a_an

d_frame <- matrix(nrow = (length(a_all)+length(an_all)+length(mic_all)),ncol = 2) 
counter <- 1
for(entry in a_all){
  d_frame[counter,1] <- entry
  d_frame[counter,2] <- 'Aerob'
  counter <- counter +1
}
for(entry in an_all){
  d_frame[counter,1] <- entry
  d_frame[counter,2] <- 'Anaerob'
  counter <- counter +1
}
for(entry in mic_all){
  d_frame[counter,1] <- entry
  d_frame[counter,2] <- 'Microaerob'
  counter <- counter +1
}


write.csv(d_frame,file = 'all.csv')

#####################################################################################################
# RNA - Group
### SETTINGS
# set (full path) the directory containing the counts (without dash at the end)
dir = "C:/Users/seide/Dropbox/Studium/Master/Masterarbeit/c_glutamicum_atcc_13032/DiffExp/counts/RNA" 
nr_datasets = 6 # number of datasets


results_dir = "C:/Users/seide/Dropbox/Studium/Master/Masterarbeit/c_glutamicum_atcc_13032/DiffExp/results/RNA"
setwd(results_dir)

# set names/labels
ntimepoints <- c("3h","5h","7h","9h","11h","13h") # timepoints
timepoints <- list(list())
timepoints[[1]]

ngroups <- c("aerob","microaerob","anaerob") # groups
groups <- list(list())
groups[[1]] <- c(1,2) # aerob
groups[[2]] <- c(3,4,5) # microaerob
groups[[3]] <- c(6) # anaerob
names(groups) <- ngroups

### DATA PREPARATION
fnames <- list.files(path=dir,full.names=T,include.dirs = FALSE) # retrieve abs path of count files
fnames <- fnames[1:3]
data <- lapply(fnames,read.table, sep="\t", header=TRUE) # read count files

types <- basename(fnames)
types <- c("asRNAs","mRNAs","ncRNAs")
types <- rep(types,sapply(data,nrow))


# merge count files (list of lists) into one
data <- do.call("rbind",data)

mylength <- data.frame(data$Geneid,data$Length)

counts <- data.frame(data$JL.3h.2.Stuttgart_CGATGT_L001.JL.3h.sorted.bam,data$JL.5h.4.Stuttgart_TGACCA_L001.JL.5h.sorted.bam,data$JL.7h.7.Stuttgart_CAGATC_L001.JL.7h.sorted.bam,data$JL.9h.16.Stuttgart_CCGTCC_L001.JL.9h.sorted.bam,data$JL.11h.18.Stuttgart_ATTCCT_L001.JL.11.sorted.bam,data$JL.13h.27.Stuttgart_GTCCGC_L001.JL.13.sorted.bam)
rownames(counts) <-data$Geneid

mybiotypes <- data.frame(types)
rownames(mybiotypes) <- data$Geneid
col_names <- c('3h','5h','7h','9h','11h','13h')
colnames(counts) <- col_names

myfactors <- data.frame(Organism = c('c_glutamicum','c_glutamicum','c_glutamicum','c_glutamicum','c_glutamicum','c_glutamicum'),Group = c('aerob-Reference', 'aerob', 'microaerob', 'microaerob', 'microaerob', 'anaerob'))

mydata <- readData(data = counts, length = mylength, factors = myfactors, biotype = mybiotypes)

# Quality Controll [!!!!!!]
mycountsbio = dat(mydata, factor = NULL, type = "countsbio")
mysaturation <- dat(mydata, k = 0, ndepth = 7, type = "saturation")
mylengthbias <- dat(mydata, factor = "Group", type = "lengthbias")
myPCA <- dat(mydata, type = "PCA")

png('explorative.png',width=10, height = 27, units = 'cm',  res= 300)
nf <- layout(matrix(c(1,2,3,4),ncol=1), widths=c(4,4,4,4), heights=c(4,4,4,4), TRUE)


explo.plot(mysaturation, toplot = 1, samples = 1:6, yleftlim = NULL, yrightlim = NULL)

explo.plot(mylengthbias, samples = NULL, toplot = "global")

explo.plot(myPCA, factor = "Group")
explo.plot(mycountsbio, toplot = 1, samples = NULL, plottype = "boxplot")
dev.off()

#Steps included in noiseq function
# Normalisation with TMM

#myTMM <- tmm(assayData(mydata)$exprs,long = mylength, lc = 1)

# filtering of Genes with low counts -> noise reduction
#myfilt <- filtered.data(counts, factor = myfactors$Group, norm = FALSE, depth = NULL, method = 1, cv.cutoff = 100, cpm = 1,p.adj = 'fdr')
prob_threshold <- 0.92
# Calculating results aerob vs aerob
myresults_a_a <- noiseq(mydata, factor = 'Group', conditions = c('aerob-Reference', 'aerob'), k = NULL, norm = 'tmm',pnr = 0.2, v = 0.02, lc = 1, replicates = 'no')
myresults_a_a.deg <- degenes(myresults_a_a, q = prob_threshold, M = NULL)
myresults_a_a.deg1 <- degenes(myresults_a_a, q = prob_threshold, M = 'up')
myresults_a_a.deg2 <- degenes(myresults_a_a, q = prob_threshold, M = 'down')


# Calculating results aerob vs anaerob
myresults_a_an <- noiseq(mydata, factor = 'Group', conditions = c('aerob-Reference', 'anaerob'), k = NULL, norm = 'tmm',pnr = 0.2, v = 0.02, lc = 1, replicates = 'no')
myresults_a_an.deg <- degenes(myresults_a_an, q = prob_threshold, M = NULL)
myresults_a_an.deg1 <- degenes(myresults_a_an, q = prob_threshold, M = 'up')
myresults_a_an.deg2 <- degenes(myresults_a_an, q = prob_threshold, M = 'down')


# Calculating results aerob vs microaerob
myresults_a_mic <- noiseq(mydata, factor = 'Group', conditions = c('aerob-Reference', 'microaerob'), k = NULL, norm = 'tmm',pnr = 0.2, v = 0.02, lc = 1, replicates = 'no')
myresults_a_mic.deg <- degenes(myresults_a_mic, q = prob_threshold, M = NULL)
myresults_a_mic.deg1 <- degenes(myresults_a_mic, q = prob_threshold, M = 'up')
myresults_a_mic.deg2 <- degenes(myresults_a_mic, q = prob_threshold, M = 'down')


# Calculating results anaerob vs microaerob
#myresults_an_mic <- noiseq(mydata, factor = 'Group', conditions = c('microaerob', 'anaerob'), k = NULL, norm = 'tmm',pnr = 0.2, v = 0.02, lc = 1, replicates = 'no')
#myresults_an_mic.deg <- degenes(myresults_an_mic, q = prob_threshold, M = NULL)
#myresults_an_mic.deg1 <- degenes(myresults_an_mic, q = prob_threshold, M = 'up')
#myresults_an_mic.deg2 <- degenes(myresults_an_mic, q = prob_threshold, M = 'down')

# Plotting of the MD DIstribution with highlighting of the DE genes
png('MD.png',width=10, height = 27, units = 'cm',  res= 300)
nf <- layout(matrix(c(1,2,3,4),ncol=1), widths=c(4,4,4,4), heights=c(4,4,4,4), TRUE)
DE.plot(myresults_a_a, q = prob_threshold, graphic = 'MD')
DE.plot(myresults_a_mic, q = prob_threshold, graphic = 'MD')
DE.plot(myresults_a_an, q = prob_threshold, graphic = 'MD')
#DE.plot(myresults_an_mic, q = prob_threshold, graphic = 'MD')
dev.off()
png('Expression.png',width=15, height = 27, units = 'cm',  res= 300)
nf<-layout(matrix(c(1,2,3,4),ncol=1), widths=c(4,4,4,4), heights=c(4,4,4,4), TRUE)
DE.plot(myresults_a_a, q = prob_threshold, graphic = "expr", log.scale = TRUE)
DE.plot(myresults_a_an, q = prob_threshold, graphic = "expr", log.scale = TRUE)
DE.plot(myresults_a_mic, q = prob_threshold, graphic = "expr", log.scale = TRUE)
#DE.plot(myresults_an_mic, q = prob_threshold, graphic = "expr", log.scale = TRUE)
dev.off()
png('Distribution.png',width=10, height = 27, units = 'cm',  res= 300)
nf<-layout(matrix(c(1,2,3,4),ncol=1), widths=c(4,4,4,4), heights=c(4,4,4,4), TRUE)
DE.plot(myresults_a_a, chromosomes = NULL, q = prob_threshold, graphic = "distr")
DE.plot(myresults_a_an, chromosomes = NULL, q = prob_threshold, graphic = "distr")
DE.plot(myresults_a_mic, chromosomes = NULL, q = prob_threshold, graphic = "distr")
#DE.plot(myresults_an_mic, chromosomes = NULL, q = prob_threshold, graphic = "distr")
dev.off()
# extract the feature names of the differential expressed features between the conditions
dif_exp_a_a <- rownames(myresults_a_a.deg)
dif_exp_a_an <- rownames(myresults_a_an.deg)
dif_exp_a_mic <- rownames(myresults_a_mic.deg)
#dif_exp_an_mic <- rownames(myresults_an_mic.deg)

# extract the feature names of the differential expressed features between the conditions
dif_exp_a_a_up <- rownames(myresults_a_a.deg1)
dif_exp_a_an_up <- rownames(myresults_a_an.deg1)
dif_exp_a_mic_up <- rownames(myresults_a_mic.deg1)
#dif_exp_an_mic_up <- rownames(myresults_an_mic.deg1)

# extract the feature names of the differential expressed features between the conditions
dif_exp_a_a_down <- rownames(myresults_a_a.deg2)
dif_exp_a_an_down <- rownames(myresults_a_an.deg2)
dif_exp_a_mic_down <- rownames(myresults_a_mic.deg2)
#dif_exp_an_mic_down <- rownames(myresults_an_mic.deg2)



dif_exp_all <- unique(c(dif_exp_a_a,dif_exp_a_an,dif_exp_a_mic))
dif_exp_all_up <- unique(c(dif_exp_a_a_up,dif_exp_a_an_up,dif_exp_a_mic_up))
dif_exp_all_down <- unique(c(dif_exp_a_a_down,dif_exp_a_an_down,dif_exp_a_mic_down))



length_df <- as.numeric(max(length(dif_exp_a_a_down),length(c(dif_exp_a_mic_down)),length(c(dif_exp_a_an_down))))
a_down <- dif_exp_a_a_down
mic_down <- c(dif_exp_a_mic_down)
an_down <- c(dif_exp_a_an_down)


d_frame <- matrix(nrow = (length(a_down)+length(an_down)+length(mic_down)),ncol = 2) 
counter <- 1
for(entry in a_down){
  d_frame[counter,1] <- entry
  d_frame[counter,2] <- 'Aerob'
  counter <- counter +1
}
for(entry in an_down){
  d_frame[counter,1] <- entry
  d_frame[counter,2] <- 'Anaerob'
  counter <- counter +1
}
for(entry in mic_down){
  d_frame[counter,1] <- entry
  d_frame[counter,2] <- 'Microaerob'
  counter <- counter +1
}


write.csv(d_frame,file = 'downregulated.csv')



# Upregulated 
a_up <- dif_exp_a_a_up
mic_up <- dif_exp_a_mic_up
an_up <- dif_exp_a_an_up

d_frame <- matrix(nrow = (length(a_up)+length(an_up)+length(mic_up)),ncol = 2) 
counter <- 1
for(entry in a_up){
  d_frame[counter,1] <- entry
  d_frame[counter,2] <- 'Aerob'
  counter <- counter +1
}
for(entry in an_up){
  d_frame[counter,1] <- entry
  d_frame[counter,2] <- 'Anaerob'
  counter <- counter +1
}
for(entry in mic_up){
  d_frame[counter,1] <- entry
  d_frame[counter,2] <- 'Microaerob'
  counter <- counter +1
}


write.csv(d_frame,file = 'upregulated.csv')

# Upregulated & Downregulated
a_all <- dif_exp_a_a
mic_all <- dif_exp_a_mic
an_all <- dif_exp_a_an

d_frame <- matrix(nrow = (length(a_all)+length(an_all)+length(mic_all)),ncol = 2) 
counter <- 1
for(entry in a_all){
  d_frame[counter,1] <- entry
  d_frame[counter,2] <- 'Aerob'
  counter <- counter +1
}
for(entry in an_all){
  d_frame[counter,1] <- entry
  d_frame[counter,2] <- 'Anaerob'
  counter <- counter +1
}
for(entry in mic_all){
  d_frame[counter,1] <- entry
  d_frame[counter,2] <- 'Microaerob'
  counter <- counter +1
}


write.csv(d_frame,file = 'all.csv')

##
library('NOISeq')
library("gplots")
library("networkD3")
#####################################################################################################
# sRNA - TP
### SETTINGS
# set (full path) the directory containing the counts (without dash at the end)
dir = "C:/Users/seide/Dropbox/Studium/Master/Masterarbeit/c_glutamicum_atcc_13032/DiffExp/counts/sRNA" 
nr_datasets = 6 # number of datasets


results_dir = "C:/Users/seide/Dropbox/Studium/Master/Masterarbeit/c_glutamicum_atcc_13032/DiffExp/results/sRNA"
setwd(results_dir)


### DATA PREPARATION
fnames <- list.files(path=dir,full.names=T,include.dirs = FALSE) # retrieve abs path of count files
fnames <- fnames[1:3]
data <- lapply(fnames,read.table, sep="\t", header=TRUE) # read count files

types <- basename(fnames)
types <- c("asRNAs","mRNAs","ncRNAs")
types <- rep(types,sapply(data,nrow))


# merge count files (list of lists) into one
data <- do.call("rbind",data)

mylength <- data.frame(data$Geneid,data$Length)

counts <- data.frame(data$X.data.jannik.data.rnaseq.sRNA_focus.JL1_sRNA_3h_AGTCAA_L001_R1_001.mapped.sorted.bam,data$X.data.jannik.data.rnaseq.sRNA_focus.JL2_sRNA_5h_AGTTCC_L001_R1_001.mapped.sorted.bam,data$X.data.jannik.data.rnaseq.sRNA_focus.JL3_sRNA_7h_ATGTCA_L001_R1_001.mapped.sorted.bam,data$X.data.jannik.data.rnaseq.sRNA_focus.JL4_sRNA_9h_CCGTCC_L001_R1_001.mapped.sorted.bam,data$X.data.jannik.data.rnaseq.sRNA_focus.JL5_sRNA_11h_GTAGAG_L001_R1_001.mapped.sorted.bam,data$X.data.jannik.data.rnaseq.sRNA_focus.JL6_sRNA_13h_GTCCGC_L001_R1_001.mapped.sorted.bam)
rownames(counts) <-data$Geneid
col_names <- c('3h','5h','7h','9h','11h','13h')
colnames(counts) <- col_names

mybiotypes <- data.frame(types)
rownames(mybiotypes) <- data$Geneid

myfactors <- data.frame(Organism = c('c_glutamicum','c_glutamicum','c_glutamicum','c_glutamicum','c_glutamicum','c_glutamicum'),Group = c('aerob-Reference', 'aerob', 'microaerob_1', 'microaerob_2', 'microaerob_3', 'anaerob'))

mydata <- readData(data = counts, length = mylength, factors = myfactors, biotype = mybiotypes)


# Quality Controll
mycountsbio = dat(mydata, factor = NULL, type = "countsbio")
mycountsbio = dat(mydata, factor = NULL, type = "countsbio")
mysaturation <- dat(mydata, k = 0, ndepth = 7, type = "saturation")
mylengthbias <- dat(mydata, factor = "Group", type = "lengthbias")
myPCA <- dat(mydata, type = "PCA")

png('explorative_tp.png',width=10, height = 27, units = 'cm',  res= 300)
nf <- layout(matrix(c(1,2,3,4),ncol=1), widths=c(4,4,4,4), heights=c(4,4,4,4), TRUE)
explo.plot(mysaturation, toplot = 1, samples = 1:6, yleftlim = NULL, yrightlim = NULL)
explo.plot(mylengthbias, samples = NULL, toplot = "global")
explo.plot(myPCA, factor = "Group")
explo.plot(mycountsbio, toplot = 1, samples = NULL, plottype = "boxplot")
dev.off()

#Steps included in noiseq function
# Normalisation with TMM

#myTMM <- tmm(assayData(mydata)$exprs,long = mylength, lc = 1)

# filtering of Genes with low counts -> noise reduction
#myfilt <- filtered.data(counts, factor = myfactors$Group, norm = FALSE, depth = NULL, method = 1, cv.cutoff = 100, cpm = 1,p.adj = 'fdr')
prob_threshold <- 0.92
# Calculating results aerob vs aerob
myresults_a_a <- noiseq(mydata, factor = 'Group', conditions = c('aerob-Reference', 'aerob'), k = NULL, norm = 'tmm',pnr = 0.2, v = 0.02, lc = 1, replicates = 'no')
myresults_a_a.deg <- degenes(myresults_a_a, q = prob_threshold, M = NULL)
myresults_a_a.deg1 <- degenes(myresults_a_a, q = prob_threshold, M = 'up')
myresults_a_a.deg2 <- degenes(myresults_a_a, q = prob_threshold, M = 'down')


# Calculating results aerob vs anaerob
myresults_a_an <- noiseq(mydata, factor = 'Group', conditions = c('aerob-Reference', 'anaerob'), k = NULL, norm = 'tmm',pnr = 0.2, v = 0.02, lc = 1, replicates = 'no')
myresults_a_an.deg <- degenes(myresults_a_an, q = prob_threshold, M = NULL)
myresults_a_an.deg1 <- degenes(myresults_a_an, q = prob_threshold, M = 'up')
myresults_a_an.deg2 <- degenes(myresults_a_an, q = prob_threshold, M = 'down')


# Calculating results aerob vs microaerob_1
myresults_a_mic_1 <- noiseq(mydata, factor = 'Group', conditions = c('aerob-Reference', 'microaerob_1'), k = NULL, norm = 'tmm',pnr = 0.2, v = 0.02, lc = 1, replicates = 'no')
myresults_a_mic_1.deg <- degenes(myresults_a_mic_1, q = prob_threshold, M = NULL)
myresults_a_mic_1.deg1 <- degenes(myresults_a_mic_1, q = prob_threshold, M = 'up')
myresults_a_mic_1.deg2 <- degenes(myresults_a_mic_1, q = prob_threshold, M = 'down')

# Calculating results aerob vs microaerob_2
myresults_a_mic_2 <- noiseq(mydata, factor = 'Group', conditions = c('aerob-Reference', 'microaerob_2'), k = NULL, norm = 'tmm',pnr = 0.2, v = 0.02, lc = 1, replicates = 'no')
myresults_a_mic_2.deg <- degenes(myresults_a_mic_2, q = prob_threshold, M = NULL)
myresults_a_mic_2.deg1 <- degenes(myresults_a_mic_2, q = prob_threshold, M = 'up')
myresults_a_mic_2.deg2 <- degenes(myresults_a_mic_2, q = prob_threshold, M = 'down')

# Calculating results aerob vs microaerob_3
myresults_a_mic_3 <- noiseq(mydata, factor = 'Group', conditions = c('aerob-Reference', 'microaerob_3'), k = NULL, norm = 'tmm',pnr = 0.2, v = 0.02, lc = 1, replicates = 'no')
myresults_a_mic_3.deg <- degenes(myresults_a_mic_3, q = prob_threshold, M = NULL)
myresults_a_mic_3.deg1 <- degenes(myresults_a_mic_3, q = prob_threshold, M = 'up')
myresults_a_mic_3.deg2 <- degenes(myresults_a_mic_3, q = prob_threshold, M = 'down')
# Calculating results anaerob vs microaerob
#myresults_an_mic <- noiseq(mydata, factor = 'Group', conditions = c('microaerob', 'anaerob'), k = NULL, norm = 'tmm',pnr = 0.2, v = 0.02, lc = 1, replicates = 'no')
#myresults_an_mic.deg <- degenes(myresults_an_mic, q = prob_threshold, M = NULL)
#myresults_an_mic.deg1 <- degenes(myresults_an_mic, q = prob_threshold, M = 'up')
#myresults_an_mic.deg2 <- degenes(myresults_an_mic, q = prob_threshold, M = 'down')

png('MD_tp.png',width=20, height = 27, units = 'cm',  res= 300)
nf <- layout(matrix(c(1,2,3,4,5,6),ncol=2), widths=c(4,4,4,4,4,0), heights=c(4,4,4,4,4,0), TRUE)
DE.plot(myresults_a_a, q = prob_threshold, graphic = 'MD')
DE.plot(myresults_a_mic_1, q = prob_threshold, graphic = 'MD')
DE.plot(myresults_a_mic_2, q = prob_threshold, graphic = 'MD')
DE.plot(myresults_a_mic_3, q = prob_threshold, graphic = 'MD')
DE.plot(myresults_a_an, q = prob_threshold, graphic = 'MD')
#DE.plot(myresults_an_mic, q = prob_threshold, graphic = 'MD')
dev.off()
png('Expression_tp.png',width=20, height = 27, units = 'cm',  res= 300)
nf <- layout(matrix(c(1,2,3,4,5,6),ncol=2), widths=c(4,4,4,4,4,0), heights=c(4,4,4,4,4,0), TRUE)
DE.plot(myresults_a_a, q = prob_threshold, graphic = "expr", log.scale = TRUE)
DE.plot(myresults_a_an, q = prob_threshold, graphic = "expr", log.scale = TRUE)
DE.plot(myresults_a_mic_1, q = prob_threshold, graphic = "expr", log.scale = TRUE)
DE.plot(myresults_a_mic_2, q = prob_threshold, graphic = "expr", log.scale = TRUE)
DE.plot(myresults_a_mic_3, q = prob_threshold, graphic = "expr", log.scale = TRUE)
#DE.plot(myresults_an_mic, q = prob_threshold, graphic = "expr", log.scale = TRUE)
dev.off()
png('Distribution_tp.png',width=20, height = 27, units = 'cm',  res= 300)
nf <- layout(matrix(c(1,2,3,4,5,6),ncol=2), widths=c(4,4,4,4,4,0), heights=c(4,4,4,4,4,0), TRUE)
DE.plot(myresults_a_a, chromosomes = NULL, q = prob_threshold, graphic = "distr")
DE.plot(myresults_a_an, chromosomes = NULL, q = prob_threshold, graphic = "distr")
DE.plot(myresults_a_mic_1, chromosomes = NULL, q = prob_threshold, graphic = "distr")
DE.plot(myresults_a_mic_2, chromosomes = NULL, q = prob_threshold, graphic = "distr")
DE.plot(myresults_a_mic_3, chromosomes = NULL, q = prob_threshold, graphic = "distr")
#DE.plot(myresults_an_mic, chromosomes = NULL, q = prob_threshold, graphic = "distr")
dev.off()
# extract the feature names of the differential expressed features between the conditions
dif_exp_a_a <- rownames(myresults_a_a.deg)
dif_exp_a_an <- rownames(myresults_a_an.deg)
dif_exp_a_mic_1 <- rownames(myresults_a_mic_1.deg)
dif_exp_a_mic_2 <- rownames(myresults_a_mic_2.deg)
dif_exp_a_mic_3 <- rownames(myresults_a_mic_3.deg)
#dif_exp_an_mic <- rownames(myresults_an_mic.deg)

# extract the feature names of the differential expressed features between the conditions
dif_exp_a_a_up <- rownames(myresults_a_a.deg1)
dif_exp_a_an_up <- rownames(myresults_a_an.deg1)
dif_exp_a_mic_up_1 <- rownames(myresults_a_mic_1.deg1)
dif_exp_a_mic_up_2 <- rownames(myresults_a_mic_2.deg1)
dif_exp_a_mic_up_3 <- rownames(myresults_a_mic_3.deg1)
#dif_exp_an_mic_up <- rownames(myresults_an_mic.deg1)

# extract the feature names of the differential expressed features between the conditions
dif_exp_a_a_down <- rownames(myresults_a_a.deg2)
dif_exp_a_an_down <- rownames(myresults_a_an.deg2)
dif_exp_a_mic_down_1 <- rownames(myresults_a_mic_1.deg2)
dif_exp_a_mic_down_2 <- rownames(myresults_a_mic_2.deg2)
dif_exp_a_mic_down_3 <- rownames(myresults_a_mic_3.deg2)
#dif_exp_an_mic_down <- rownames(myresults_an_mic.deg2)



dif_exp_all <- unique(c(dif_exp_a_a,dif_exp_a_an,dif_exp_a_mic_1,dif_exp_a_mic_2,dif_exp_a_mic_3))
dif_exp_all_up <- unique(c(dif_exp_a_a_up,dif_exp_a_an_up,dif_exp_a_mic_up_1,dif_exp_a_mic_up_2,dif_exp_a_mic_up_3))
dif_exp_all_down <- unique(c(dif_exp_a_a_down,dif_exp_a_an_down,dif_exp_a_mic_down_1, dif_exp_a_mic_down_2,dif_exp_a_mic_down_3))



length_df <- as.numeric(max(length(dif_exp_a_a_down),length(c(dif_exp_a_mic_down_1)),length(c(dif_exp_a_mic_down_2)),length(c(dif_exp_a_mic_down_3)),length(c(dif_exp_a_an_down))))
a_down <- dif_exp_a_a_down
mic_down_1 <- c(dif_exp_a_mic_down_1)
mic_down_2 <- c(dif_exp_a_mic_down_2)
mic_down_3 <- c(dif_exp_a_mic_down_3)
an_down <- c(dif_exp_a_an_down)


d_frame <- matrix(nrow = (length(a_down)+length(an_down)+length(mic_down_1)+length(mic_down_2)+length(mic_down_3)),ncol = 2) 
counter <- 1
for(entry in a_down){
  d_frame[counter,1] <- entry
  d_frame[counter,2] <- 'Aerob'
  counter <- counter +1
}
for(entry in an_down){
  d_frame[counter,1] <- entry
  d_frame[counter,2] <- 'Anaerob'
  counter <- counter +1
}
for(entry in mic_down_1){
  d_frame[counter,1] <- entry
  d_frame[counter,2] <- 'Microaerob_1'
  counter <- counter +1
}
for(entry in mic_down_2){
  d_frame[counter,1] <- entry
  d_frame[counter,2] <- 'Microaerob_2'
  counter <- counter +1
}
for(entry in mic_down_3){
  d_frame[counter,1] <- entry
  d_frame[counter,2] <- 'Microaerob_3'
  counter <- counter +1
}
write.csv(d_frame,file = 'downregulated_tp.csv')



# Upregulated 
a_up <- dif_exp_a_a_up
mic_up_1 <- dif_exp_a_mic_up_1
mic_up_2 <- dif_exp_a_mic_up_2
mic_up_3 <- dif_exp_a_mic_up_3
an_up <- dif_exp_a_an_up

d_frame <- matrix(nrow = (length(a_up)+length(an_up)+length(mic_up_1)+length(mic_up_2)+length(mic_up_3)),ncol = 2) 
counter <- 1
for(entry in a_up){
  d_frame[counter,1] <- entry
  d_frame[counter,2] <- 'Aerob'
  counter <- counter +1
}
for(entry in an_up){
  d_frame[counter,1] <- entry
  d_frame[counter,2] <- 'Anaerob'
  counter <- counter +1
}
for(entry in mic_up_1){
  d_frame[counter,1] <- entry
  d_frame[counter,2] <- 'Microaerob_1'
  counter <- counter +1
}
for(entry in mic_up_2){
  d_frame[counter,1] <- entry
  d_frame[counter,2] <- 'Microaerob_2'
  counter <- counter +1
}
for(entry in mic_up_3){
  d_frame[counter,1] <- entry
  d_frame[counter,2] <- 'Microaerob_3'
  counter <- counter +1
}

write.csv(d_frame,file = 'upregulated_tp.csv')


# Upregulated & Downregulated
a_all <- dif_exp_a_a
mic_all_1 <- dif_exp_a_mic_1
mic_all_2 <- dif_exp_a_mic_2
mic_all_3 <- dif_exp_a_mic_3
an_all <- dif_exp_a_an

d_frame <- matrix(nrow = (length(a_all)+length(an_all)+length(mic_all_1)+length(mic_all_2)+length(mic_all_3)),ncol = 2) 
counter <- 1
for(entry in a_all){
  d_frame[counter,1] <- entry
  d_frame[counter,2] <- 'Aerob'
  counter <- counter +1
}
for(entry in an_all){
  d_frame[counter,1] <- entry
  d_frame[counter,2] <- 'Anaerob'
  counter <- counter +1
}
for(entry in mic_all_1){
  d_frame[counter,1] <- entry
  d_frame[counter,2] <- 'Microaerob_1'
  counter <- counter +1
}
for(entry in mic_all_2){
  d_frame[counter,1] <- entry
  d_frame[counter,2] <- 'Microaerob_2'
  counter <- counter +1
}
for(entry in mic_all_3){
  d_frame[counter,1] <- entry
  d_frame[counter,2] <- 'Microaerob_3'
  counter <- counter +1
}



write.csv(d_frame,file = 'all_tp.csv')

names <- unique(d_frame[,1])
data_frame <- data.frame(matrix(ncol = 5, nrow = length(names)), row.names = names)
col_names <- c('aerob','microaerob1','microaerob2','microaerob3','anaerob')
colnames(data_frame) <- col_names
# construct chord diagram from dataframe
for (entry in a_all){
  data_frame[entry,1]<- 1
  
}
for (entry in mic_all_1){
  data_frame[entry,2]<- 1
  
}
for (entry in mic_all_2){
  data_frame[entry,3]<- 1
  
}
for (entry in mic_all_3){
  data_frame[entry,4]<- 1
  
}
for (entry in an_all){
  data_frame[entry,5]<- 1
}

adj_list <- data.frame(source=character(),target=character())
for (i in 1:length(data_frame$aerob)){
  
}


#####################################################################################################
# RNA - TP
### SETTINGS
# set (full path) the directory containing the counts (without dash at the end)
dir = "C:/Users/seide/Dropbox/Studium/Master/Masterarbeit/c_glutamicum_atcc_13032/DiffExp/counts/RNA" 
nr_datasets = 6 # number of datasets


results_dir = "C:/Users/seide/Dropbox/Studium/Master/Masterarbeit/c_glutamicum_atcc_13032/DiffExp/results/RNA"
setwd(results_dir)

# set names/labels
ntimepoints <- c("3h","5h","7h","9h","11h","13h") # timepoints
timepoints <- list(list())
timepoints[[1]]

ngroups <- c("aerob","microaerob","anaerob") # groups
groups <- list(list())
groups[[1]] <- c(1,2) # aerob
groups[[2]] <- c(3,4,5) # microaerob
groups[[3]] <- c(6) # anaerob
names(groups) <- ngroups

### DATA PREPARATION
fnames <- list.files(path=dir,full.names=T,include.dirs = FALSE) # retrieve abs path of count files
fnames <- fnames[1:3]
data <- lapply(fnames,read.table, sep="\t", header=TRUE) # read count files

types <- basename(fnames)
types <- c("asRNAs","mRNAs","ncRNAs")
types <- rep(types,sapply(data,nrow))


# merge count files (list of lists) into one
data <- do.call("rbind",data)

mylength <- data.frame(data$Geneid,data$Length)

counts <- data.frame(data$JL.3h.2.Stuttgart_CGATGT_L001.JL.3h.sorted.bam,data$JL.5h.4.Stuttgart_TGACCA_L001.JL.5h.sorted.bam,data$JL.7h.7.Stuttgart_CAGATC_L001.JL.7h.sorted.bam,data$JL.9h.16.Stuttgart_CCGTCC_L001.JL.9h.sorted.bam,data$JL.11h.18.Stuttgart_ATTCCT_L001.JL.11.sorted.bam,data$JL.13h.27.Stuttgart_GTCCGC_L001.JL.13.sorted.bam)
rownames(counts) <-data$Geneid

mybiotypes <- data.frame(types)
rownames(mybiotypes) <- data$Geneid
col_names <- c('3h','5h','7h','9h','11h','13h')
colnames(counts) <- col_names

myfactors <- data.frame(Organism = c('c_glutamicum','c_glutamicum','c_glutamicum','c_glutamicum','c_glutamicum','c_glutamicum'),Group = c('aerob-Reference', 'aerob', 'microaerob_1', 'microaerob_2', 'microaerob_3', 'anaerob'))

mydata <- readData(data = counts, length = mylength, factors = myfactors, biotype = mybiotypes)

# Quality Controll [!!!!!!]
mycountsbio = dat(mydata, factor = NULL, type = "countsbio")
mysaturation <- dat(mydata, k = 0, ndepth = 7, type = "saturation")
mylengthbias <- dat(mydata, factor = "Group", type = "lengthbias")
myPCA <- dat(mydata, type = "PCA")

png('explorative_tp.png',width=10, height = 27, units = 'cm',  res= 300)
nf <- layout(matrix(c(1,2,3,4),ncol=1), widths=c(4,4,4,4), heights=c(4,4,4,4), TRUE)


explo.plot(mysaturation, toplot = 1, samples = 1:6, yleftlim = NULL, yrightlim = NULL)

explo.plot(mylengthbias, samples = NULL, toplot = "global")

explo.plot(myPCA, factor = "Group")
explo.plot(mycountsbio, toplot = 1, samples = NULL, plottype = "boxplot")
dev.off()

#Steps included in noiseq function
# Normalisation with TMM

#myTMM <- tmm(assayData(mydata)$exprs,long = mylength, lc = 1)

# filtering of Genes with low counts -> noise reduction
#myfilt <- filtered.data(counts, factor = myfactors$Group, norm = FALSE, depth = NULL, method = 1, cv.cutoff = 100, cpm = 1,p.adj = 'fdr')
prob_threshold <- 0.92
# Calculating results aerob vs aerob
myresults_a_a <- noiseq(mydata, factor = 'Group', conditions = c('aerob-Reference', 'aerob'), k = NULL, norm = 'tmm',pnr = 0.2, v = 0.02, lc = 1, replicates = 'no')
myresults_a_a.deg <- degenes(myresults_a_a, q = prob_threshold, M = NULL)
myresults_a_a.deg1 <- degenes(myresults_a_a, q = prob_threshold, M = 'up')
myresults_a_a.deg2 <- degenes(myresults_a_a, q = prob_threshold, M = 'down')


# Calculating results aerob vs anaerob
myresults_a_an <- noiseq(mydata, factor = 'Group', conditions = c('aerob-Reference', 'anaerob'), k = NULL, norm = 'tmm',pnr = 0.2, v = 0.02, lc = 1, replicates = 'no')
myresults_a_an.deg <- degenes(myresults_a_an, q = prob_threshold, M = NULL)
myresults_a_an.deg1 <- degenes(myresults_a_an, q = prob_threshold, M = 'up')
myresults_a_an.deg2 <- degenes(myresults_a_an, q = prob_threshold, M = 'down')


# Calculating results aerob vs microaerob_1
myresults_a_mic_1 <- noiseq(mydata, factor = 'Group', conditions = c('aerob-Reference', 'microaerob_1'), k = NULL, norm = 'tmm',pnr = 0.2, v = 0.02, lc = 1, replicates = 'no')
myresults_a_mic_1.deg <- degenes(myresults_a_mic_1, q = prob_threshold, M = NULL)
myresults_a_mic_1.deg1 <- degenes(myresults_a_mic_1, q = prob_threshold, M = 'up')
myresults_a_mic_1.deg2 <- degenes(myresults_a_mic_1, q = prob_threshold, M = 'down')

# Calculating results aerob vs microaerob_2
myresults_a_mic_2 <- noiseq(mydata, factor = 'Group', conditions = c('aerob-Reference', 'microaerob_2'), k = NULL, norm = 'tmm',pnr = 0.2, v = 0.02, lc = 1, replicates = 'no')
myresults_a_mic_2.deg <- degenes(myresults_a_mic_2, q = prob_threshold, M = NULL)
myresults_a_mic_2.deg1 <- degenes(myresults_a_mic_2, q = prob_threshold, M = 'up')
myresults_a_mic_2.deg2 <- degenes(myresults_a_mic_2, q = prob_threshold, M = 'down')

# Calculating results aerob vs microaerob_3
myresults_a_mic_3 <- noiseq(mydata, factor = 'Group', conditions = c('aerob-Reference', 'microaerob_3'), k = NULL, norm = 'tmm',pnr = 0.2, v = 0.02, lc = 1, replicates = 'no')
myresults_a_mic_3.deg <- degenes(myresults_a_mic_3, q = prob_threshold, M = NULL)
myresults_a_mic_3.deg1 <- degenes(myresults_a_mic_3, q = prob_threshold, M = 'up')
myresults_a_mic_3.deg2 <- degenes(myresults_a_mic_3, q = prob_threshold, M = 'down')
# Calculating results anaerob vs microaerob
#myresults_an_mic <- noiseq(mydata, factor = 'Group', conditions = c('microaerob', 'anaerob'), k = NULL, norm = 'tmm',pnr = 0.2, v = 0.02, lc = 1, replicates = 'no')
#myresults_an_mic.deg <- degenes(myresults_an_mic, q = prob_threshold, M = NULL)
#myresults_an_mic.deg1 <- degenes(myresults_an_mic, q = prob_threshold, M = 'up')
#myresults_an_mic.deg2 <- degenes(myresults_an_mic, q = prob_threshold, M = 'down')

# Plotting of the MD DIstribution with highlighting of the DE genes
png('MD_tp.png',width=20, height = 27, units = 'cm',  res= 300)
nf <- layout(matrix(c(1,2,3,4,5,6),ncol=2), widths=c(4,4,4,4,4,0), heights=c(4,4,4,4,4,0), TRUE)
DE.plot(myresults_a_a, q = prob_threshold, graphic = 'MD')
DE.plot(myresults_a_mic_1, q = prob_threshold, graphic = 'MD')
DE.plot(myresults_a_mic_2, q = prob_threshold, graphic = 'MD')
DE.plot(myresults_a_mic_3, q = prob_threshold, graphic = 'MD')
DE.plot(myresults_a_an, q = prob_threshold, graphic = 'MD')
#DE.plot(myresults_an_mic, q = prob_threshold, graphic = 'MD')
dev.off()
png('Expression_tp.png',width=20, height = 27, units = 'cm',  res= 300)
nf <- layout(matrix(c(1,2,3,4,5,6),ncol=2), widths=c(4,4,4,4,4,0), heights=c(4,4,4,4,4,0), TRUE)
DE.plot(myresults_a_a, q = prob_threshold, graphic = "expr", log.scale = TRUE)
DE.plot(myresults_a_an, q = prob_threshold, graphic = "expr", log.scale = TRUE)
DE.plot(myresults_a_mic_1, q = prob_threshold, graphic = "expr", log.scale = TRUE)
DE.plot(myresults_a_mic_2, q = prob_threshold, graphic = "expr", log.scale = TRUE)
DE.plot(myresults_a_mic_3, q = prob_threshold, graphic = "expr", log.scale = TRUE)
#DE.plot(myresults_an_mic, q = prob_threshold, graphic = "expr", log.scale = TRUE)
dev.off()
png('Distribution_tp.png',width=20, height = 27, units = 'cm',  res= 300)
nf <- layout(matrix(c(1,2,3,4,5,6),ncol=2), widths=c(4,4,4,4,4,0), heights=c(4,4,4,4,4,0), TRUE)
DE.plot(myresults_a_a, chromosomes = NULL, q = prob_threshold, graphic = "distr")
DE.plot(myresults_a_an, chromosomes = NULL, q = prob_threshold, graphic = "distr")
DE.plot(myresults_a_mic_1, chromosomes = NULL, q = prob_threshold, graphic = "distr")
DE.plot(myresults_a_mic_2, chromosomes = NULL, q = prob_threshold, graphic = "distr")
DE.plot(myresults_a_mic_3, chromosomes = NULL, q = prob_threshold, graphic = "distr")
#DE.plot(myresults_an_mic, chromosomes = NULL, q = prob_threshold, graphic = "distr")
dev.off()
# extract the feature names of the differential expressed features between the conditions
dif_exp_a_a <- rownames(myresults_a_a.deg)
dif_exp_a_an <- rownames(myresults_a_an.deg)
dif_exp_a_mic_1 <- rownames(myresults_a_mic_1.deg)
dif_exp_a_mic_2 <- rownames(myresults_a_mic_2.deg)
dif_exp_a_mic_3 <- rownames(myresults_a_mic_3.deg)
#dif_exp_an_mic <- rownames(myresults_an_mic.deg)

# extract the feature names of the differential expressed features between the conditions
dif_exp_a_a_up <- rownames(myresults_a_a.deg1)
dif_exp_a_an_up <- rownames(myresults_a_an.deg1)
dif_exp_a_mic_up_1 <- rownames(myresults_a_mic_1.deg1)
dif_exp_a_mic_up_2 <- rownames(myresults_a_mic_2.deg1)
dif_exp_a_mic_up_3 <- rownames(myresults_a_mic_3.deg1)
#dif_exp_an_mic_up <- rownames(myresults_an_mic.deg1)

# extract the feature names of the differential expressed features between the conditions
dif_exp_a_a_down <- rownames(myresults_a_a.deg2)
dif_exp_a_an_down <- rownames(myresults_a_an.deg2)
dif_exp_a_mic_down_1 <- rownames(myresults_a_mic_1.deg2)
dif_exp_a_mic_down_2 <- rownames(myresults_a_mic_2.deg2)
dif_exp_a_mic_down_3 <- rownames(myresults_a_mic_3.deg2)
#dif_exp_an_mic_down <- rownames(myresults_an_mic.deg2)



dif_exp_all <- unique(c(dif_exp_a_a,dif_exp_a_an,dif_exp_a_mic_1,dif_exp_a_mic_2,dif_exp_a_mic_3))
dif_exp_all_up <- unique(c(dif_exp_a_a_up,dif_exp_a_an_up,dif_exp_a_mic_up_1,dif_exp_a_mic_up_2,dif_exp_a_mic_up_3))
dif_exp_all_down <- unique(c(dif_exp_a_a_down,dif_exp_a_an_down,dif_exp_a_mic_down_1, dif_exp_a_mic_down_2,dif_exp_a_mic_down_3))



length_df <- as.numeric(max(length(dif_exp_a_a_down),length(c(dif_exp_a_mic_down_1)),length(c(dif_exp_a_mic_down_2)),length(c(dif_exp_a_mic_down_3)),length(c(dif_exp_a_an_down))))
a_down <- dif_exp_a_a_down
mic_down_1 <- c(dif_exp_a_mic_down_1)
mic_down_2 <- c(dif_exp_a_mic_down_2)
mic_down_3 <- c(dif_exp_a_mic_down_3)
an_down <- c(dif_exp_a_an_down)


d_frame <- matrix(nrow = (length(a_down)+length(an_down)+length(mic_down_1)+length(mic_down_2)+length(mic_down_3)),ncol = 2) 
counter <- 1
for(entry in a_down){
  d_frame[counter,1] <- entry
  d_frame[counter,2] <- 'Aerob'
  counter <- counter +1
}
for(entry in an_down){
  d_frame[counter,1] <- entry
  d_frame[counter,2] <- 'Anaerob'
  counter <- counter +1
}
for(entry in mic_down_1){
  d_frame[counter,1] <- entry
  d_frame[counter,2] <- 'Microaerob_1'
  counter <- counter +1
}
for(entry in mic_down_2){
  d_frame[counter,1] <- entry
  d_frame[counter,2] <- 'Microaerob_2'
  counter <- counter +1
}
for(entry in mic_down_3){
  d_frame[counter,1] <- entry
  d_frame[counter,2] <- 'Microaerob_3'
  counter <- counter +1
}
write.csv(d_frame,file = 'downregulated_tp.csv')
names <- unique(d_frame[,1])
data_frame <- data.frame(matrix(ncol = 5, nrow = length(names)), row.names = names)
col_names <- c('aerob','microaerob1','microaerob2','microaerob3','anaerob')
colnames(data_frame) <- col_names
for (entry in a_down){
  data_frame[entry,1]<- 1
  
}
for (entry in mic_down_1){
  data_frame[entry,2]<- 1
  
}
for (entry in mic_down_2){
  data_frame[entry,3]<- 1
  
}
for (entry in mic_down_3){
  data_frame[entry,4]<- 1
  
}
for (entry in an_down){
  data_frame[entry,5]<- 1
  
}


data_frame[is.na(data_frame)]<-0


# Upregulated 
a_up <- dif_exp_a_a_up
mic_up_1 <- dif_exp_a_mic_up_1
mic_up_2 <- dif_exp_a_mic_up_2
mic_up_3 <- dif_exp_a_mic_up_3
an_up <- dif_exp_a_an_up

d_frame <- matrix(nrow = (length(a_up)+length(an_up)+length(mic_up_1)+length(mic_up_2)+length(mic_up_3)),ncol = 2) 
counter <- 1
for(entry in a_up){
  d_frame[counter,1] <- entry
  d_frame[counter,2] <- 'Aerob'
  counter <- counter +1
}
for(entry in an_up){
  d_frame[counter,1] <- entry
  d_frame[counter,2] <- 'Anaerob'
  counter <- counter +1
}
for(entry in mic_up_1){
  d_frame[counter,1] <- entry
  d_frame[counter,2] <- 'Microaerob_1'
  counter <- counter +1
}
for(entry in mic_up_2){
  d_frame[counter,1] <- entry
  d_frame[counter,2] <- 'Microaerob_2'
  counter <- counter +1
}
for(entry in mic_up_3){
  d_frame[counter,1] <- entry
  d_frame[counter,2] <- 'Microaerob_3'
  counter <- counter +1
}

write.csv(d_frame,file = 'upregulated_tp.csv')

names <- unique(d_frame[,1])
data_frame <- data.frame(matrix(ncol = 5, nrow = length(names)), row.names = names)
col_names <- c('aerob','microaerob1','microaerob2','microaerob3','anaerob')
colnames(data_frame) <- col_names

for (entry in a_up){
  data_frame[entry,1]<- 1
  
}
for (entry in mic_up_1){
  data_frame[entry,2]<- 1
  
}
for (entry in mic_up_2){
  data_frame[entry,3]<- 1
  
}
for (entry in mic_up_3){
  data_frame[entry,4]<- 1
  
}
for (entry in an_up){
  data_frame[entry,5]<- 1
  
}

data_frame[is.na(data_frame)]<-0

# Upregulated & Downregulated
a_all <- dif_exp_a_a
mic_all_1 <- dif_exp_a_mic_1
mic_all_2 <- dif_exp_a_mic_2
mic_all_3 <- dif_exp_a_mic_3
an_all <- dif_exp_a_an

d_frame <- matrix(nrow = (length(a_all)+length(an_all)+length(mic_all_1)+length(mic_all_2)+length(mic_all_3)),ncol = 2) 
counter <- 1
for(entry in a_all){
  d_frame[counter,1] <- entry
  d_frame[counter,2] <- 'Aerob'
  counter <- counter +1
}
for(entry in an_all){
  d_frame[counter,1] <- entry
  d_frame[counter,2] <- 'Anaerob'
  counter <- counter +1
}
for(entry in mic_all_1){
  d_frame[counter,1] <- entry
  d_frame[counter,2] <- 'Microaerob_1'
  counter <- counter +1
}
for(entry in mic_all_2){
  d_frame[counter,1] <- entry
  d_frame[counter,2] <- 'Microaerob_2'
  counter <- counter +1
}
for(entry in mic_all_3){
  d_frame[counter,1] <- entry
  d_frame[counter,2] <- 'Microaerob_3'
  counter <- counter +1
}



write.csv(d_frame,file = 'all_tp.csv')

names <- unique(d_frame[,1])
data_frame <- data.frame(matrix(ncol = 5, nrow = length(names)), row.names = names)
col_names <- c('aerob','microaerob1','microaerob2','microaerob3','anaerob')
colnames(data_frame) <- col_names
# construct chord diagram from dataframe
for (entry in a_all){
  data_frame[entry,1]<- 1
  
}
for (entry in mic_all_1){
  data_frame[entry,2]<- 1
  
}
for (entry in mic_all_2){
  data_frame[entry,3]<- 1
  
}
for (entry in mic_all_3){
  data_frame[entry,4]<- 1
  
}
for (entry in an_all){
  data_frame[entry,5]<- 1
  
}

data_frame[is.na(data_frame)]<-0

pydir = "C:/Users/seide/Dropbox/Studium/Master/Masterarbeit/c_glutamicum_atcc_13032/DiffExp/"
setwd(pydir)

py_run_file("find_name_start_stop_strand.py")
py_run_file("extract_as_nc_RNA.py")
py_run_file("find_as.py")


# RNA - TP Network
library(readr)
library(networkD3)
library(plyr)
library(stringi)
downregulated_tp_feature <- read_delim("C:/Users/seide/Dropbox/Studium/Master/Masterarbeit/c_glutamicum_atcc_13032/DiffExp/results/RNA/downregulated_tp_feature.txt", 
                                       "\t", escape_double = FALSE, col_names = FALSE, 
                                       trim_ws = TRUE)
down <- rep('down',length(downregulated_tp_feature$X1))
downregulated_tp_feature$direction <- down
upregulated_tp_feature <- read_delim("C:/Users/seide/Dropbox/Studium/Master/Masterarbeit/c_glutamicum_atcc_13032/DiffExp/results/RNA/upregulated_tp_feature.txt", 
                                     "\t", escape_double = FALSE, col_names = FALSE, 
                                     trim_ws = TRUE)
up <- rep('up', length(upregulated_tp_feature$X1))
upregulated_tp_feature$direction <- up

names <- unique(c(downregulated_tp_feature$X1,upregulated_tp_feature$X1))
nodes <- data.frame(name=c("Aerob","Microaerob_1", "Microaerob_2","Microaerob_3","Anaerob"),group='timepoint',regulation='',stringsAsFactors=F)
counter = 1
for (i in downregulated_tp_feature$X1){
  if (!(i %in% nodes$name)&& substr(i,1,1) == 'c'){
    len_regulation <- table(unlist(all_tp$X1))[[i]]
    regulation <- character(len_regulation[1])
    sub_frame <- all_tp[all_tp$X1==i,]
    for (j in 1:length(sub_frame$X1)){
      regulation[j] <- paste(sub_frame[j,'X6'],sub_frame[j,'direction'],sep = ':')

    }
    regulation <- paste(regulation,sep = '; ',collapse = '; ')
    df <- data.frame(name = i, group=downregulated_tp_feature$X5[counter],regulation=regulation,stringsAsFactors=F)
    nodes <- rbind(nodes,df)
  }
  counter <- counter +1
  
}
counter = 1
for (i in upregulated_tp_feature$X1){
  if (!(i %in% nodes$name) && substr(i,1,1) == 'c'){
    len_regulation <- table(unlist(all_tp$X1))[[i]]
    regulation <- character(len_regulation[1])
    sub_frame <- all_tp[all_tp$X1==i,]
    for (j in 1:length(sub_frame$X1)){
      regulation[j] <- paste(sub_frame[j,'X6'],sub_frame[j,'direction'],sep = ':')
    }
    regulation <- paste(regulation,sep = '; ',collapse = '; ')
    df <- data.frame(name = i, group=upregulated_tp_feature$X5[counter],regulation=regulation,stringsAsFactors=F)
    nodes <- rbind(nodes,df)
  }
  counter <- counter +1
  
}

link <- data.frame(source=numeric(),target=numeric())

all_tp <- rbind(downregulated_tp_feature,upregulated_tp_feature)
all_tp <- all_tp[order(all_tp$X1),]
counter = 0
sum(all_tp == 'cgb_00165')
for (i in nodes$name){
  if (!(i=='Aerob') && !(i=='Anaerob') && !(i=='Microaerob_1')&& !(i=='Microaerob_2')&& !(i=='Microaerob_3')){
    pos <- match(i,all_tp$X1)
    last <- pos + sum(all_tp == i)-1
    line = all_tp[pos,]
    cond <- line$X6
    pos_cond <- match(cond,nodes$name)-1
    df <- data.frame(source=counter,target=pos_cond)
    
    if (!(pos == last)){
    for (j in seq(pos+1,last,by=1)){
      line = all_tp[j,]
      if (j == last)
      {
        print('next')
      }
      if (is.na(line$X1)){
        break
      }
      cond <- line$X6
      pos_cond <- match(cond,nodes$name)-1
      df <- rbind(df,data.frame(source=counter,target=pos_cond))
    }
    }
    link <- rbind(link,df)
    
  }
  
  counter <- counter +1
  print(counter)
  print('next1------------------------')
}
new_nodes <- nodes
nodes <- new_nodes
for (i in seq(6,length(nodes$name))){
  print(i)
  new_name <- paste('id: ',nodes$name[i],sep = '')
  print(new_name)
  new_reg <- paste('regulation: ',nodes$regulation[i],sep = '')
  print(new_reg)
  new_name <- paste(new_name, new_reg, sep = ';')
  print(new_name)
  print(nodes[i,1])
  nodes$name[i]<- new_name
  print(nodes$name[i])
}

MyClickScript <- 'alert("You clicked " + d.name + " which is in group " +  d.group);'

forceNetwork(Links = link,Nodes = nodes, NodeID = 'name', Group = 'group',opacity = 0.9, zoom = TRUE, clickAction = MyClickScript)

# sRNA - TP Network

downregulated_tp_feature <- read_delim("C:/Users/seide/Dropbox/Studium/Master/Masterarbeit/c_glutamicum_atcc_13032/DiffExp/results/sRNA/downregulated_tp_feature.txt", 
                                       "\t", escape_double = FALSE, col_names = FALSE, 
                                       trim_ws = TRUE)
down <- rep('down',length(downregulated_tp_feature$X1))
downregulated_tp_feature$direction <- down
upregulated_tp_feature <- read_delim("C:/Users/seide/Dropbox/Studium/Master/Masterarbeit/c_glutamicum_atcc_13032/DiffExp/results/sRNA/upregulated_tp_feature.txt", 
                                     "\t", escape_double = FALSE, col_names = FALSE, 
                                     trim_ws = TRUE)
up <- rep('up', length(upregulated_tp_feature$X1))
upregulated_tp_feature$direction <- up

all_tp <- rbind(downregulated_tp_feature,upregulated_tp_feature)
all_tp <- all_tp[order(all_tp$X1),]

names <- unique(c(downregulated_tp_feature$X1,upregulated_tp_feature$X1))
nodes <- data.frame(name=c("Aerob","Microaerob_1", "Microaerob_2","Microaerob_3","Anaerob"),group='timepoint',regulation='',stringsAsFactors=F)
counter = 1
for (i in downregulated_tp_feature$X1){
  if (!(i %in% nodes$name)&& substr(i,1,1) == 'c'){
    len_regulation <- table(unlist(all_tp$X1))[[i]]
    regulation <- character(len_regulation[1])
    sub_frame <- all_tp[all_tp$X1==i,]
    for (j in 1:length(sub_frame$X1)){
      regulation[j] <- paste(sub_frame[j,'X6'],sub_frame[j,'direction'],sep = ':')
      
    }
    regulation <- paste(regulation,sep = '; ',collapse = '; ')
    df <- data.frame(name = i, group=downregulated_tp_feature$X5[counter],regulation=regulation,stringsAsFactors=F)
    nodes <- rbind(nodes,df)
  }
  counter <- counter +1
  
}
counter = 1
for (i in upregulated_tp_feature$X1){
  if (!(i %in% nodes$name) && substr(i,1,1) == 'c'){
    len_regulation <- table(unlist(all_tp$X1))[[i]]
    regulation <- character(len_regulation[1])
    sub_frame <- all_tp[all_tp$X1==i,]
    for (j in 1:length(sub_frame$X1)){
      regulation[j] <- paste(sub_frame[j,'X6'],sub_frame[j,'direction'],sep = ':')
    }
    regulation <- paste(regulation,sep = '; ',collapse = '; ')
    df <- data.frame(name = i, group=upregulated_tp_feature$X5[counter],regulation=regulation,stringsAsFactors=F)
    nodes <- rbind(nodes,df)
  }
  counter <- counter +1
  
}

link <- data.frame(source=numeric(),target=numeric())


counter = 0
sum(all_tp == 'cgb_00165')
for (i in nodes$name){
  if (!(i=='Aerob') && !(i=='Anaerob') && !(i=='Microaerob_1')&& !(i=='Microaerob_2')&& !(i=='Microaerob_3')){
    pos <- match(i,all_tp$X1)
    last <- pos + sum(all_tp == i)-1
    line = all_tp[pos,]
    cond <- line$X6
    pos_cond <- match(cond,nodes$name)-1
    df <- data.frame(source=counter,target=pos_cond)
    
    if (!(pos == last)){
      for (j in seq(pos+1,last,by=1)){
        line = all_tp[j,]
        if (j == last)
        {
          print('next')
        }
        if (is.na(line$X1)){
          break
        }
        cond <- line$X6
        pos_cond <- match(cond,nodes$name)-1
        df <- rbind(df,data.frame(source=counter,target=pos_cond))
      }
    }
    link <- rbind(link,df)
    
  }
  
  counter <- counter +1
  print(counter)
  print('next1------------------------')
}
new_nodes <- nodes
nodes <- new_nodes
for (i in seq(6,length(nodes$name))){
  print(i)
  new_name <- paste('id: ',nodes$name[i],sep = '')
  print(new_name)
  new_reg <- paste('regulation: ',nodes$regulation[i],sep = '')
  print(new_reg)
  new_name <- paste(new_name, new_reg, sep = ';')
  print(new_name)
  print(nodes[i,1])
  nodes$name[i]<- new_name
  print(nodes$name[i])
}

MyClickScript <- 'alert("You clicked " + d.name + " which is in group " +  d.group);'

forceNetwork(Links = link,Nodes = nodes, NodeID = 'name', Group = 'group',opacity = 0.9, zoom = TRUE, clickAction = MyClickScript)
#---------------------------------------------------------------------------------------------------
# Interaktionsnetzwerk
# RNA-TP
downregulated_tp_feature <- read_delim("C:/Users/seide/Dropbox/Studium/Master/Masterarbeit/c_glutamicum_atcc_13032/DiffExp/results/RNA/downregulated_tp_feature.txt", 
                                       "\t", escape_double = FALSE, col_names = FALSE, 
                                       trim_ws = TRUE)
down <- rep('down',length(downregulated_tp_feature$X1))
downregulated_tp_feature$direction <- down
upregulated_tp_feature <- read_delim("C:/Users/seide/Dropbox/Studium/Master/Masterarbeit/c_glutamicum_atcc_13032/DiffExp/results/RNA/upregulated_tp_feature.txt", 
                                     "\t", escape_double = FALSE, col_names = FALSE, 
                                     trim_ws = TRUE)
up <- rep('up', length(upregulated_tp_feature$X1))
upregulated_tp_feature$direction <- up

down_as <- read_delim("C:/Users/seide/Dropbox/Studium/Master/Masterarbeit/c_glutamicum_atcc_13032/DiffExp/results/RNA/downregulated_tp_feature_as.txt", 
                      "\t", escape_double = FALSE, col_names = FALSE, 
                      trim_ws = TRUE)
up_as <- read_delim("C:/Users/seide/Dropbox/Studium/Master/Masterarbeit/c_glutamicum_atcc_13032/DiffExp/results/RNA/upregulated_tp_feature_as.txt", 
                    "\t", escape_double = FALSE, col_names = FALSE, 
                    trim_ws = TRUE)
all_tp <- rbind(downregulated_tp_feature,upregulated_tp_feature)
all_tp <- all_tp[order(all_tp$X1),]
all_as <- rbind(down_as,up_as)
all_as <- subset(all_as,select = -c(X6))
all_as <- unique(all_as)
for (i in 1:length(all_as$X1)){
  str_as <- all_as$X5[i]
  str_as <- unlist(strsplit(str_as,';',fixed = TRUE))[2]
  str_as <- unlist(strsplit(str_as,'=', fixed = TRUE))[2]
  all_as$X5[i] <- str_as
}




names <- unique(c(downregulated_tp_feature$X1,upregulated_tp_feature$X1))
nodes <- data.frame(name=c("Aerob","Microaerob_1", "Microaerob_2","Microaerob_3","Anaerob",'Missing'),group=c('timepoint','timepoint','timepoint','timepoint','timepoint','missing'),regulation='',stringsAsFactors=F)

counter = 1
for (i in downregulated_tp_feature$X1){
  if (!(i %in% nodes$name)){
    len_regulation <- table(unlist(all_tp$X1))[[i]]
    regulation <- character(len_regulation[1])
    sub_frame <- all_tp[all_tp$X1==i,]
    for (j in 1:length(sub_frame$X1)){
      regulation[j] <- paste(sub_frame[j,'X6'],sub_frame[j,'direction'],sep = ':')
      
    }
    regulation <- paste(regulation,sep = '; ',collapse = '; ')
    df <- data.frame(name = i, group=downregulated_tp_feature$X5[counter],regulation=regulation,stringsAsFactors=F)
    nodes <- rbind(nodes,df)
  }
  counter <- counter +1
  
}
counter = 1
for (i in upregulated_tp_feature$X1){
  if (!(i %in% nodes$name) ){
    len_regulation <- table(unlist(all_tp$X1))[[i]]
    regulation <- character(len_regulation[1])
    sub_frame <- all_tp[all_tp$X1==i,]
    for (j in 1:length(sub_frame$X1)){
      regulation[j] <- paste(sub_frame[j,'X6'],sub_frame[j,'direction'],sep = ':')
    }
    regulation <- paste(regulation,sep = '; ',collapse = '; ')
    df <- data.frame(name = i, group=upregulated_tp_feature$X5[counter],regulation=regulation,stringsAsFactors=F)
    nodes <- rbind(nodes,df)
  }
  counter <- counter +1
  
}

for (i in all_as$X5){
  if (!(i %in% nodes$name)){
    df <- data.frame(name = i, group='missing',regulation='',stringsAsFactors=F)
    nodes <- rbind(nodes,df)
  }
}

link <- data.frame(source=numeric(),target=numeric())


counter = 0
sum(all_tp == 'cgb_00165')
for (i in nodes$name){
  df_defined <- FALSE
  if (!(i=='Aerob') && !(i=='Anaerob') && !(i=='Microaerob_1')&& !(i=='Microaerob_2')&& !(i=='Microaerob_3') && !(i=='Missing')){
    pos <- match(i,all_tp$X1)
    last <- pos + sum(all_tp == i)-1
    line = all_tp[pos,]
    cond <- line$X6
    pos_cond <- match(cond,nodes$name)-1
    if (!(is.na(pos))){
      df <- data.frame(source=counter,target=pos_cond)
      df_defined <- TRUE
    
      if (!(is.na(pos)) && !(pos == last)){
        for (j in seq(pos+1,last,by=1)){
          line = all_tp[j,]
          if (j == last)
          {
            print('next')
          }
          if (is.na(line$X1)){
            break
          }
          cond <- line$X6
          pos_cond <- match(cond,nodes$name)-1
          df <- rbind(df,data.frame(source=counter,target=pos_cond))
        }}}
    if (i %in% all_as$X1){
      
      targ <- match(i, all_as$X1)
      targ <- all_as$X5[targ]
      targ <- match(targ,nodes$name)-1
      print('TARGET:')
      print(targ)
      if (df_defined){
        df <- rbind(df, data.frame(source=counter, target = targ))}
      else{
        df <-  data.frame(source=counter, target = targ)
        df_defined <- TRUE
      }
    }
    
    if (i %in% nodes$name && nodes$group[match(i,nodes$name)] == 'missing' && !(i == 'Missing')){
      targe <- 6-1
      if (df_defined){
        df <- rbind(df,data.frame(source = counter, target = targe))
      }
      else{
        df <- data.frame(source = counter, target = targe)
        df_defined <- TRUE
      }
    }
    if (df_defined){
      link <- rbind(link,df)
    }
  }
  
  counter <- counter +1
  print(counter)
  print('next1------------------------')
}
new_nodes <- nodes
nodes <- new_nodes
for (i in seq(6,length(nodes$name))){
  print(i)
  new_name <- paste('id: ',nodes$name[i],sep = '')
  print(new_name)
  new_reg <- paste('regulation: ',nodes$regulation[i],sep = '')
  print(new_reg)
  new_name <- paste(new_name, new_reg, sep = ';')
  print(new_name)
  print(nodes[i,1])
  nodes$name[i]<- new_name
  print(nodes$name[i])
}

MyClickScript <- 'alert("You clicked " + d.name + " which is in group " +  d.group);'

forceNetwork(Links = link,Nodes = nodes, NodeID = 'name', Group = 'group',opacity = 0.9, zoom = TRUE, clickAction = MyClickScript)

# sRNA - TP Network

downregulated_tp_feature <- read_delim("C:/Users/seide/Dropbox/Studium/Master/Masterarbeit/c_glutamicum_atcc_13032/DiffExp/results/sRNA/downregulated_tp_feature.txt", 
                                       "\t", escape_double = FALSE, col_names = FALSE, 
                                       trim_ws = TRUE)
down <- rep('down',length(downregulated_tp_feature$X1))
downregulated_tp_feature$direction <- down
upregulated_tp_feature <- read_delim("C:/Users/seide/Dropbox/Studium/Master/Masterarbeit/c_glutamicum_atcc_13032/DiffExp/results/sRNA/upregulated_tp_feature.txt", 
                                     "\t", escape_double = FALSE, col_names = FALSE, 
                                     trim_ws = TRUE)
up <- rep('up', length(upregulated_tp_feature$X1))
upregulated_tp_feature$direction <- up

down_as <- read_delim("C:/Users/seide/Dropbox/Studium/Master/Masterarbeit/c_glutamicum_atcc_13032/DiffExp/results/sRNA/downregulated_tp_feature_as.txt", 
                      "\t", escape_double = FALSE, col_names = FALSE, 
                      trim_ws = TRUE)
up_as <- read_delim("C:/Users/seide/Dropbox/Studium/Master/Masterarbeit/c_glutamicum_atcc_13032/DiffExp/results/sRNA/upregulated_tp_feature_as.txt", 
                    "\t", escape_double = FALSE, col_names = FALSE, 
                    trim_ws = TRUE)
all_tp <- rbind(downregulated_tp_feature,upregulated_tp_feature)
all_tp <- all_tp[order(all_tp$X1),]
all_as <- rbind(down_as,up_as)
all_as <- subset(all_as,select = -c(X6))
all_as <- unique(all_as)
for (i in 1:length(all_as$X1)){
  str_as <- all_as$X5[i]
  str_as <- unlist(strsplit(str_as,';',fixed = TRUE))[2]
  str_as <- unlist(strsplit(str_as,'=', fixed = TRUE))[2]
  all_as$X5[i] <- str_as
}




names <- unique(c(downregulated_tp_feature$X1,upregulated_tp_feature$X1))
nodes <- data.frame(name=c("Aerob","Microaerob_1", "Microaerob_2","Microaerob_3","Anaerob",'Missing'),group=c('timepoint','timepoint','timepoint','timepoint','timepoint','missing'),regulation='',stringsAsFactors=F)

counter = 1
for (i in downregulated_tp_feature$X1){
  if (!(i %in% nodes$name)){
    len_regulation <- table(unlist(all_tp$X1))[[i]]
    regulation <- character(len_regulation[1])
    sub_frame <- all_tp[all_tp$X1==i,]
    for (j in 1:length(sub_frame$X1)){
      regulation[j] <- paste(sub_frame[j,'X6'],sub_frame[j,'direction'],sep = ':')
      
    }
    regulation <- paste(regulation,sep = '; ',collapse = '; ')
    df <- data.frame(name = i, group=downregulated_tp_feature$X5[counter],regulation=regulation,stringsAsFactors=F)
    nodes <- rbind(nodes,df)
  }
  counter <- counter +1
  
}
counter = 1
for (i in upregulated_tp_feature$X1){
  if (!(i %in% nodes$name) ){
    len_regulation <- table(unlist(all_tp$X1))[[i]]
    regulation <- character(len_regulation[1])
    sub_frame <- all_tp[all_tp$X1==i,]
    for (j in 1:length(sub_frame$X1)){
      regulation[j] <- paste(sub_frame[j,'X6'],sub_frame[j,'direction'],sep = ':')
    }
    regulation <- paste(regulation,sep = '; ',collapse = '; ')
    df <- data.frame(name = i, group=upregulated_tp_feature$X5[counter],regulation=regulation,stringsAsFactors=F)
    nodes <- rbind(nodes,df)
  }
  counter <- counter +1
  
}

for (i in all_as$X5){
  if (!(i %in% nodes$name)){
    df <- data.frame(name = i, group='missing',regulation='',stringsAsFactors=F)
    nodes <- rbind(nodes,df)
  }
}

link <- data.frame(source=numeric(),target=numeric())


counter = 0
sum(all_tp == 'cgb_00165')
for (i in nodes$name){
  df_defined <- FALSE
  if (!(i=='Aerob') && !(i=='Anaerob') && !(i=='Microaerob_1')&& !(i=='Microaerob_2')&& !(i=='Microaerob_3') && !(i=='Missing')){
    pos <- match(i,all_tp$X1)
    last <- pos + sum(all_tp == i)-1
    line = all_tp[pos,]
    cond <- line$X6
    pos_cond <- match(cond,nodes$name)-1
    if (!(is.na(pos))){
      df <- data.frame(source=counter,target=pos_cond)
      df_defined <- TRUE
      
      if (!(is.na(pos)) && !(pos == last)){
        for (j in seq(pos+1,last,by=1)){
          line = all_tp[j,]
          if (j == last)
          {
            print('next')
          }
          if (is.na(line$X1)){
            break
          }
          cond <- line$X6
          pos_cond <- match(cond,nodes$name)-1
          df <- rbind(df,data.frame(source=counter,target=pos_cond))
        }}}
    if (i %in% all_as$X1){
      
      targ <- match(i, all_as$X1)
      targ <- all_as$X5[targ]
      targ <- match(targ,nodes$name)-1
      print('TARGET:')
      print(targ)
      if (df_defined){
        df <- rbind(df, data.frame(source=counter, target = targ))}
      else{
        df <-  data.frame(source=counter, target = targ)
        df_defined <- TRUE
      }
    }
    
    if (i %in% nodes$name && nodes$group[match(i,nodes$name)] == 'missing' && !(i == 'Missing')){
      targe <- 6-1
      if (df_defined){
        df <- rbind(df,data.frame(source = counter, target = targe))
      }
      else{
        df <- data.frame(source = counter, target = targe)
        df_defined <- TRUE
      }
    }
    if (df_defined){
      link <- rbind(link,df)
    }
  }
  
  counter <- counter +1
  print(counter)
  print('next1------------------------')
}


new_nodes <- nodes
nodes <- new_nodes
for (i in seq(6,length(nodes$name))){
  print(i)
  new_name <- paste('id: ',nodes$name[i],sep = '')
  print(new_name)
  new_reg <- paste('regulation: ',nodes$regulation[i],sep = '')
  print(new_reg)
  new_name <- paste(new_name, new_reg, sep = ';')
  print(new_name)
  print(nodes[i,1])
  nodes$name[i]<- new_name
  print(nodes$name[i])
}

MyClickScript <- 'alert("You clicked " + d.name + " which is in group " +  d.group);'

forceNetwork(Links = link,Nodes = nodes, NodeID = 'name', Group = 'group',opacity = 0.9, zoom = TRUE, clickAction = MyClickScript)






