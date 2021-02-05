#!/bin/bash
# script for mapping the trimmed reads to the index and sorting and converting to bam

for DIR in ls -d /data/jannik/data/rnaseq/sRNA_focus/*;
do	
	echo 'STARTED'
        
	name1="${DIR}/output.fastq"
	sam="$DIR/mapped.sam"
	echo $sam
	sam1="$DIR/mapped.bam"
	bam="$DIR/mapped.sorted.bam"
	bai="$DIR/mapped.sorted.bam.bai"
        bowtie2 -x /data/jannik/data/rnaseq/bowtie2index/c_glutamicum1tcc13032 -U ${name1} -S ${sam}
	samtools view -S -b ${sam}>${sam1}
	samtools sort ${sam1} > ${bam}
	samtools index ${bam} > ${bai}
	echo 'HEUREKA'
    

        
done
echo 'FINISHED'
