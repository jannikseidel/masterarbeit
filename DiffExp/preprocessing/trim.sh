#!/bin/bash
# script for trimming the single end reads from different conditions

for DIR in ls -d /data/jannik/data/rnaseq/sRNA_focus/*;
do 
	
	array=()
	counter=0
	for file in "$DIR"*;
	do
		array[$counter]=$file
		let "counter=counter+1"
	done

	file=$(ls ${array[0]})

	echo "${array[0]}/${file[0]}"
	java -jar /data/jannik/data/rnaseq/tools/Trimmomatic-0.38/trimmomatic-0.38.jar SE -phred33 "${array[0]}/${file[0]}" "${array[0]}/output.fastq" ILLUMINACLIP:/data/c_glutamicum/JulianRNA-Sequencing/procreads/adapters.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
done
