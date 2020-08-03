#!/bin/bash

# run it within the experiment
# pass all MID numbers as arguments

MIDS="$@"

echo "MIDS: $MIDS"

echo "making stats files.."

for i in $MIDS
do
    fastx_quality_stats -i *MID${i}_*_R1_*001.fastq -o mid${i}_r1_stats.txt
    fastx_quality_stats -i *MID${i}_*_R2_*001.fastq -o mid${i}_r2_stats.txt
done

echo "making final statistics.."

for r in 1 2
do
    echo "MIDs #Reads AvgMeans AvgQ1s AvgMedians AvgQ3s" > r${r}_stats.csv
    for i in $MIDS
    do
        echo -n "MID${i} "
        tail -n+2 mid${i}_r${r}_stats.txt | awk '{means+=$6;Q1s+=$7;meds+=$8;Q3s+=$9} END {print $2" "means/NR" "Q1s/NR" "meds/NR" "Q3s/NR}'
    done >> r${r}_stats.csv
done

rm *_stats.txt
