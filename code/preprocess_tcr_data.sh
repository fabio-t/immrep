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

# echo "trimming.."

# for f in *.fastq
# do
# 	fastx_trimmer -f 1 -l 150 -i $f -o `basename $f .fastq`_trimmed.fastq
# done

# echo "making stats for trimmed files.."

# for i in $MIDS
# do
# 	fastx_quality_stats -i *MID${i}_*_R1_*001_trimmed.fastq -o mid${i}_r1_trimmed_stats.txt
# 	fastx_quality_stats -i *MID${i}_*_R2_*001_trimmed.fastq -o mid${i}_r2_trimmed_stats.txt
# done

echo "making final statistics.."

for trim in ""
do
	for r in 1 2
	do
		echo "MIDs #Reads AvgMeans AvgQ1s AvgMedians AvgQ3s" > r${r}${trim}_stats.csv
		for i in $MIDS
		do
			echo -n "MID${i} "
			tail -n+2 mid${i}_r${r}${trim}_stats.txt | awk '{means+=$6;Q1s+=$7;meds+=$8;Q3s+=$9} END {print $2" "means/NR" "Q1s/NR" "meds/NR" "Q3s/NR}'
		done >> r${r}${trim}_stats.csv
	done
done

rm *_stats.txt
