#!/bin/bash
set -e

# TODO this should be turned into functions, so scripts can compose it as they see fit

# scripts must set:
# $resdir $indices $exp_data $exp_data_dir $data_type $bc $R

# TODO add default values, provide a better output at the beginning,
# as a reminder of the various variables that can be called in experiment
# scripts

jointype=${jointype:-gene}

m1=$m
m=${m+_m${m}}

ql1=$ql
ql=${ql+_ql${ql}}

execscript=${execscript:-$exp_dir}

echo $resdir
echo $bc
echo $R
echo $tool
echo $m
echo $ql
echo $execscript
echo $exp_dir
echo ${exp_dir}_${tool}-m${m1}
echo $exp_data_dir
echo $data_type
echo $indices
echo $jointype

data_prefix="/mnt/storage/data/local/mol_med/${data_type}/${exp_data_dir}/samples"

# FIXME more checks/output here
# TODO add || exit 1 with a proper error correction after the main commands
mkdir -p ${exp_dir}_${tool}-m${m1} && cd ${exp_dir}_${tool}-m${m1}
mkdir -p stats collisions full
for i in $indices
do
    if [[ $tool == "migec" ]]
    then
      cat ${data_prefix}/${resdir}/${R}${m}${bc}${ql}/cdrfinal_${i}/S0_${R}.csv | awk -F"\t" '{print $1"\t"$3"\t"$4"\t"$5"\t"$6"\t"$2}' > mid${i}_clones.csv
      cat ${data_prefix}/${resdir}/${R}${m}${bc}${ql}/cdrfinal_${i}/S0_${R}.csv > full/MID${i}.csv
    else
      mixcr exportClones -t -o -count -nFeature CDR3 -aaFeature CDR3 -vHit -jHit -dHit -fraction -f ${data_prefix}/${resdir}/${R}${m}${bc}${ql}/mixcr_${i}/analysis.*clns mid${i}_clones.csv
      # mixcr exportClones -t -o -count -nFeature CDR3 -aaFeature CDR3 -vHit -jHit -dHit -cHit -fraction -f ${data_prefix}/${resdir}/${R}${m}${bc}${ql}/mixcr_${i}/analysis.*clns mid${i}_clones.csv
      # mixcr exportClones -t -o -count -nFeature CDR3 -aaFeature CDR3 -vGene -jGene -dGene -cGene -fraction -f ${data_prefix}/${resdir}/${R}${m}${bc}${ql}/mixcr_${i}/analysis.*clns mid${i}_clones.csv
      mixcr exportClones -t -o -vHit -jHit -dHit -cHit -p full -f ${data_prefix}/${resdir}/${R}${m}${bc}${ql}/mixcr_${i}/analysis.*clns full/MID${i}.csv
      # cp ${data_prefix}/${resdir}/${R}${m}/mixcr_${i}/*clonotypes*.txt full/MID${i}.csv
    fi
    cat mid${i}_clones.csv |  awk '{print $2}' | sort -k1,1 | uniq -c | awk '{if($1 > 1){print $0}}' | sort -k1,1 -r > collisions/mid${i}_collision.csv
done
cp ${data_prefix}/${resdir}/${R}${m}${bc}${ql}/*.csv stats/
cp ${data_prefix}/r{1,2}_stats.csv stats/
[ -f "${data_prefix}/../mid_labels.csv" ] && cp ${data_prefix}/../mid_labels.csv .

# TODO: need a script to remove duplicates in place and recalculate counts and frequencies (for cristoph's scripts)

/mnt/storage/data/code/join_mids.py --unique --type $jointype $indices

# TODO: either full or not, must be a choice
Rscript /mnt/storage/data/code/family_usage_full.R

# TODO: check if script exists? this should be automatised better, too - metadata should be pulled from somewhere
# and loaded into the R scripts
Rscript ../code/${execscript}.R
