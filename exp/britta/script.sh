#!/bin/bash

migec="migec"
safe="_safe"

# indices=`echo 16 17 19 20 {22..28} 38 {61..69}`
# exp_dir="181019_igg"
# exp_data_dir="181019_igg"

# indices=`echo {40..60}`
# exp_dir="181105_iga"
# exp_data_dir="181105_iga"

# indices=`echo {2..7} 11 {13..17} {19..21}`
# exp_dir="200623_iga"
# exp_data_dir="200623_IBD_IgA"

indices=`echo {22..40} 64`
exp_dir="200625_igg"
exp_data_dir="200625_IBD_IgG"

mkdir -p $exp_dir && cd $exp_dir

R=R12
# R=R1

mkdir -p stats collisions
for i in $indices
do
    cat ../../../local/mol_med/bcr/${exp_data_dir}/samples/${migec}/${R}/cdrfinal_${i}${safe}/S0_R12.csv | awk -F"\t" '{print $1"\t"$3"\t"$4"\t"$5"\t"$6"\t"$2}' > mid${i}_clones.csv
    cat mid${i}_clones.csv |  awk '{print $2}' | uniq -c | awk '{if($1 > 1){print $0}}' | sort -k1,1 -r > collisions/mid${i}_collision.csv
    cp ../../../local/mol_med/bcr/${exp_data_dir}/samples/${migec}/${R}/*.csv stats/
    cp ../../../local/mol_med/bcr/${exp_data_dir}/samples/r{1,2}_stats.csv stats/
done

../../../code/join_mids.py --unique --type gene $indices

Rscript ../../../code/family_usage_full.R

Rscript ../code/${exp_dir}.R
