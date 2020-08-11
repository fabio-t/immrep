#!/bin/bash

migec="migec"
safe="_safe"

indices=`echo {25..39} {5..9} {45..49} 40 21 55 56 59`
exp_dir="190218_treg"
exp_data_dir="190218_irf4flox"

indices=`echo {20..27} 37 {28..36} {38..55}`
exp_dir="200415_resident"
exp_data_dir="200415_resident"

mkdir -p $exp_dir && cd $exp_dir

R=R12
# R=R1

mkdir -p stats collisions
for i in $indices
do
    cat ../../../local/mol_med/tcr_cd3/${exp_data_dir}/samples/${migec}/${R}/cdrfinal_${i}${safe}/S0_R12.csv | awk -F"\t" '{print $1"\t"$3"\t"$4"\t"$5"\t"$6"\t"$2}' > mid${i}_clones.csv
    cat mid${i}_clones.csv |  awk '{print $2}' | sort -k1,1 | uniq -c | awk '{if($1 > 1){print $0}}' | sort -k1,1 -r > collisions/mid${i}_collision.csv
    cp ../../../local/mol_med/tcr_cd3/${exp_data_dir}/samples/${migec}/${R}/*.csv stats/
    cp ../../../local/mol_med/tcr_cd3/${exp_data_dir}/samples/r{1,2}_stats.csv stats/
done

../../../code/join_mids.py --unique --type gene $indices

Rscript ../../../code/family_usage_full.R

Rscript ../code/${exp_dir}.R
