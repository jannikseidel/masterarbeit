#!/bin/bash
# script for running MAFFT MSA on the output of GLASSgo;
# 

input_folder=/data/jannik/data/results/glassgo/*
output_folder=/data/jannik/data/results/mafft/
counter=1
for entry in $input_folder; do
	echo $counter
	file_name=${entry##*/}
	out_file=$output_folder$file_name
	in_file=$entry
	mafft --maxiterate 1000 --thread 20 $in_file > $out_file
	counter=counter+1
done
