#!/bin/bash

# TODO this should be turned into functions, so scripts can compose it as they see fit

# scripts must set:
# $migec $indices $exp_data $exp_data_dir $safe $R

echo $migec
echo $safe
echo $R
echo $exp_data
echo $exp_data_dir
echo $indices

# FIXME more checks/output here
mkdir -p $exp_dir && cd $exp_dir
mkdir -p stats collisions
for i in $indices
do
    cat ../../../local/mol_med/tcr_cd3/${exp_data_dir}/samples/${migec}/${R}/cdrfinal_${i}${safe}/S0_R12.csv | awk -F"\t" '{print $1"\t"$3"\t"$4"\t"$5"\t"$6"\t"$2}' > mid${i}_clones.csv
    cat mid${i}_clones.csv |  awk '{print $2}' | sort -k1,1 | uniq -c | awk '{if($1 > 1){print $0}}' | sort -k1,1 -r > collisions/mid${i}_collision.csv
    cp ../../../local/mol_med/tcr_cd3/${exp_data_dir}/samples/${migec}/${R}/*.csv stats/
    cp ../../../local/mol_med/tcr_cd3/${exp_data_dir}/samples/r{1,2}_stats.csv stats/
done

# TODO: this needs to be parameterised (unique or not, type gene or something else..)
../../../code/join_mids.py --unique --type gene $indices

# TODO: either full or not, must be a choice
Rscript ../../../code/family_usage_full.R

# TODO: check if script exists? this should be automatised better, too - metadata should be pulled from somewhere
# and loaded into the R scripts
Rscript ../code/${exp_dir}.R
