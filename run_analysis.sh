#!/bin/bash
# script for running the whole analysis from glassgo to mafft to CopraRNA and
# IntaRNA

# Running GLASSgo
command='python3 /data/jannik/data/scripts/glassgo/run_glassgo.py'
start=`date +%s`
$command
end=`date +%s`
time_used=$((end-start))
echo 'GLASSgo ran!'>> /data/jannik/data/results/stats/stats.txt
echo $time_used >> /data/jannik/data/results/stats/stats.txt

# Running MAFFT
command='bash /data/jannik/data/scripts/mafft/run_mafft.sh'
start=`date +%s`
$command
end=`date +%s`
time_used=$((end-start))
echo 'MAFFT ran!' >> /data/jannik/data/results/stats/stats.txt
echo $time_used >> /data/jannik/data/results/stats/stats.txt

# Filtering
command='python3 /data/jannik/data/scripts/mafft/filter.py'
start=`date +%s`
$command
end=`date +%s`
time_used=$((end-start))
echo 'Filtering ran!' >> /data/jannik/data/results/stats/stats.txt
echo $time_used >> /data/jannik/data/results/stats/stats.txt

# Examination and moving
command='python3 /data/jannik/data/scripts/mafft/examin_filtered.py'
start=`date +%s`
$command
end=`date +%s`
time_used=$((end-start))
echo 'Examination and moving finished!' >> /data/jannik/data/results/stats/stats.txt
echo $time_used >> /data/jannik/data/results/stats/stats.txt

# Running CopraRNA
command='python3 /data/jannik/data/scripts/CopraRNA/run_coprarna.py'
start=`date +%s`
$command
end=`date +%s`
time_used=$((end-start))
echo 'CopraRNA ran!' >> /data/jannik/data/results/stats/stats.txt
echo $time_used >> /data/jannik/data/results/stats/stats.txt

# Running IntaRNA
command='python3 /data/jannik/data/scripts/IntaRNA/run_intarna.py'
start=`date +%s`
$command
end=`date +%s`
time_used=$((end-start))
echo 'IntaRNA ran!' >> /data/jannik/data/results/stats/stats.txt
echo $time_used >> /data/jannik/data/results/stats/stats.txt

cat /data/jannik/data/results/stats/stats.txt
