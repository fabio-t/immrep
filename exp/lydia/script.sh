#!/bin/bash

migec="migec"
safe="_safe"

# indices=`echo {1..11} {13..17} {19..21}`
# exp_dir="170714_cultures"
# exp_data_dir="170714_cultures"

# indices=`echo {1..11} {13..17} {19..32}`
# exp_dir="171215_cultures"
# exp_data_dir="171215_cultures"

# indices=`echo {33..62}`
# exp_dir="171215_cultures_TRB"
# exp_data_dir="171215_cultures"

# indices=`echo {30..64}`
# exp_dir="181204_treg"
# exp_data_dir="181204_treg"

# indices=`echo 10 11 {13..17} 19 20 {22..24} {50..54} {62..65}`
# exp_dir="190218_irf4flox"
# exp_data_dir="190218_irf4flox"

# indices=`echo {24..46}`
# exp_dir="170912_capture"
# exp_data_dir="190320_capture"

# indices=`echo 48 50 51 52 54 55 56 58 59 57 60 62 63 61 49 53`
# exp_dir="190107_capture"
# exp_data_dir="190320_capture"

# indices=`echo 2 1 {4..11} 13 14 {25..33} {15..17} {19..24}`
# exp_dir="180830_treg"
# exp_data_dir="190715_treg"

# indices=`echo {35..44}`
# exp_dir="180830_sort"
# exp_data_dir="190715_treg"

# indices=`echo {20..55}`
# exp_dir="190220_treg"
# exp_data_dir="191104_mixed"

# indices=`echo {56..64}`
# exp_dir="190220_pools"
# exp_data_dir="191104_mixed"

# indices=`echo {27..38}`
# exp_dir="191213_itregtcon2"
# exp_data_dir="191213_mixed"

# indices=`echo {39..47}`
# exp_dir="190930_treg"
# exp_data_dir="191213_mixed"

# indices=`echo 55 {48..50} 52 53 54 51`
# exp_dir="191213_klrg1"
# exp_data_dir="191213_mixed"

# indices=`echo {20..55}`
# exp_dir="200415_resident"
# exp_data_dir="200415_resident"

# indices=`echo {1..11} {13..17} {19..47}`
# exp_dir="200615_diversity"
# exp_data_dir="200615_diversity"

## "strict" version, it requires a change in the $safe and $migec variables
# indices=`echo {35..44}`
# exp_dir="180830_sort_strict"
# exp_data_dir="190715_treg"
# migec="migec2"
# safe="_safe_long_2"

## removed chimeric and pseudogenes
indices=`echo {1..11} {13..17} {19..47}`
exp_dir="200615_diversity_2"
exp_data_dir="200615_diversity"
migec="migec2"

mkdir -p $exp_dir && cd $exp_dir

R=R12
# R=R1

mkdir -p stats collisions
for i in $indices
do
    cat ../../../local/mol_med/tcr_cd3/${exp_data_dir}/samples/${migec}/${R}/cdrfinal_${i}${safe}/S0_R12.csv | awk -F"\t" '{print $1"\t"$3"\t"$4"\t"$5"\t"$6"\t"$2}' > mid${i}_clones.csv
    cat mid${i}_clones.csv |  awk '{print $2}' | uniq -c | awk '{if($1 > 1){print $0}}' | sort -k1,1 -r > collisions/mid${i}_collision.csv
    cp ../../../local/mol_med/tcr_cd3/${exp_data_dir}/samples/${migec}/${R}/*.csv stats/
    cp ../../../local/mol_med/tcr_cd3/${exp_data_dir}/samples/r{1,2}_stats.csv stats/
done

../../../code/join_mids.py --unique --type gene $indices

Rscript ../../../code/family_usage_full.R

Rscript ../code/${exp_dir}.R
