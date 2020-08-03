#!/bin/bash

indices=`echo {35..44}`
exp_dir="180830_sort_noUMIs"
exp_data_dir="190715_treg"

mkdir -p $exp_dir && cd $exp_dir

R=R12
# R=R1

safe="_safe"
# safe=""

mkdir -p full logs
for i in $indices
do
    mixcr analyze amplicon -s mouse --starting-material rna --5-end v-primers --3-end j-primers --adapters adapters-present --receptor-type TCR --align "-t 4 ${library}" --assemble "-t 4" --export "-p full" ../../../local/mol_med/tcr_cd3/${exp_data_dir}/samples/*MID${i}_*_R{1,2}_001${trim}.fastq logs/mid${i} || exit 1
    mixcr exportClones -count -nFeature CDR3 -aaFeature CDR3 -vHit -jHit -fraction -f logs/mid${i}.clna mid${i}_clones.csv || exit 1
    cp logs/mid${i}.clonotypes.TRA.txt full/mid${i}_clones.csv || exit 1
done

../../../code/join_mids.py --unique --type gene $indices

Rscript ../../../code/family_usage.R

# Rscript ../code/${exp_dir}.R
