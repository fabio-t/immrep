#!/bin/bash

echo "MID mouse organ clones_count single_clones tot_reads aligned_reads avg_reads_per_clone reads_used_clonotypes"

for mid in "$@"
do
    stats=(`cat logs/align_mid${mid}.txt | awk -F":" 'NR>6{print $2}' | awk '{print $1}'`)

    tot_reads=${stats[0]}
    aligned_reads=${stats[1]}

    stats=(`cat logs/assemble_mid${mid}.txt | awk -F":" 'NR>6{print $2}' | awk '{print $1}'`)

    clones_count=${stats[0]}
    avg_reads_per_clone=${stats[1]}
    reads_used_clonotypes=${stats[2]}

    single_clones=`egrep "^1\s" -c mid${mid}_clones.csv`

    echo $mid X X $clones_count $single_clones $tot_reads $aligned_reads $avg_reads_per_clone $reads_used_clonotypes

done

## aligner stats
# (6 to skip)
# Total sequencing reads: 6149
# Successfully aligned reads: 5708 (92.83%)
# Alignment failed, no hits (not TCR/IG?): 86 (1.4%)
# Alignment failed because of absence of V hits: 4 (0.07%)
# Alignment failed because of absence of J hits: 339 (5.51%)
# Alignment failed because of low total score: 12 (0.2%)
# Overlapped: 5 (0.08%)
# Overlapped and aligned: 0 (0%)
# Overlapped and not aligned: 5 (0.08%)
# [TRA] chains: 5708 (100%)

## assembler stats
# (6 to skip)
# Final clonotype count: 1629
# Average number of reads per clonotype: 2.92
# Reads used in clonotypes, percent of total: 4752 (77.28%)
# Reads used in clonotypes before clustering, percent of total: 4758 (77.38%)
# Number of reads used as a core, percent of used: 4432 (93.15%)
# Mapped low quality reads, percent of used: 326 (6.85%)
# Reads clustered in PCR error correction, percent of used: 6 (0.13%)
# Reads pre-clustered due to the similar VJC-lists, percent of used: 0 (0%)
# Reads dropped due to the lack of a clone sequence: 24 (0.39%)
# Reads dropped due to low quality: 0 (0%)
# Reads dropped due to failed mapping: 926 (15.06%)
# Reads dropped with low quality clones: 0 (0%)
# Clonotypes eliminated by PCR error correction: 6
# Clonotypes dropped as low quality: 0
# Clonotypes pre-clustered due to the similar VJC-lists: 0
# [TRA] chains: 1629 (100%)
