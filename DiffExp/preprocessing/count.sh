#!/bin/bash
# script for counting the mapped reads vai featurecount

# looping over the folders
array=()
counter=1
for DIR in /data/jannik/data/rnaseq/sRNA_focus/* ;do
	folder=$DIR
	if [[ "$folder" == *"JL"* ]]; then
		for FILE in "$DIR/"*;do
			file=$FILE

			if [[ "$file" == *'.sorted.bam'* ]] && [[ "$file" != *".bai"* ]]; then
				array[$counter]=$file
				let counter+=1
				echo $counter
				echo $file
			fi
		done
	fi

done
echo ${array[6]}

featureCounts -T 5 -t gene -g ID -a /data/jannik/data/annotations/merged.gff3 -o /data/jannik/data/rnaseq/sRNA_focus/counts/counts_mRNAs.txt ${array[1]} ${array[2]} ${array[3]} ${array[4]} ${array[5]} ${array[6]}

featureCounts -T 5 -t antisense_RNA -g ID -a /data/jannik/data/annotations/merged.gff3 -o /data/jannik/data/rnaseq/sRNA_focus/counts/counts_asRNAs.txt ${array[1]} ${array[2]} ${array[3]} ${array[4]} ${array[5]} ${array[6]}

featureCounts -T 5 -t ncRNA -g ID -a /data/jannik/data/annotations/merged.gff3 -o /data/jannik/data/rnaseq/sRNA_focus/counts/counts_ncRNAs.txt ${array[1]} ${array[2]} ${array[3]} ${array[4]} ${array[5]} ${array[6]}

