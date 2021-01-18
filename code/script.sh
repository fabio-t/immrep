#!/bin/bash

# TODO this should be turned into functions, so scripts can compose it as they see fit

# scripts must set:
# $migec $indices $exp_data $exp_data_dir $data_type $safe $R

# TODO add default values, provide a better output at the beginning,
# as a reminder of the various variables that can be called in experiment
# scripts

echo $migec
echo $safe
echo $R
echo $exp_data
echo $exp_data_dir
echo $data_type
echo $indices

data_prefix="/mnt/storage/data/local/mol_med/${data_type}/${exp_data_dir}/samples"

# FIXME more checks/output here
# TODO add || exit 1 with a proper error correction after the main commands
mkdir -p $exp_dir && cd $exp_dir
mkdir -p stats collisions
for i in $indices
do
    cat ${data_prefix}/${migec}/${R}/cdrfinal_${i}${safe}/S0_R12.csv | awk -F"\t" '{print $1"\t"$3"\t"$4"\t"$5"\t"$6"\t"$2}' > mid${i}_clones.csv
    cat mid${i}_clones.csv |  awk '{print $2}' | sort -k1,1 | uniq -c | awk '{if($1 > 1){print $0}}' | sort -k1,1 -r > collisions/mid${i}_collision.csv
done
cp ${data_prefix}/${migec}/${R}/*.csv stats/
cp ${data_prefix}/r{1,2}_stats.csv stats/

# TODO: this needs to be parameterised (unique or not, type gene or something else..)
/mnt/storage/data/code/join_mids.py --unique --type gene $indices

# TODO: either full or not, must be a choice
Rscript /mnt/storage/data/code/family_usage_full.R

# TODO: check if script exists? this should be automatised better, too - metadata should be pulled from somewhere
# and loaded into the R scripts
Rscript ../code/${exp_dir}.R
