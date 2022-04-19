#!/bin/bash
set -e

# TODO this should be turned into functions, so scripts can compose it as they see fit

# scripts must set:
# $resdir $indices $exp_data $exp_data_dir $data_type $bc $R

# TODO add default values, provide a better output at the beginning,
# as a reminder of the various variables that can be called in experiment
# scripts

jointype=${jointype:-gene}

organism=${organism:-mouse}

m1=$m
m=${m+_m${m}}

ql1=$ql
ql=${ql+_ql${ql}}

execscript=${execscript:-$exp_dir}
outdir="${exp_dir}_${tool}_m${m1}${ql}"

echo $resdir
echo $bc
echo $R
echo $tool
echo $m
echo $ql
echo $execscript
echo $exp_dir
echo $outdir
echo $exp_data_dir
echo $data_type
echo $chain
echo $indices
echo $jointype
echo $downsample
echo $groupby

data_prefix="/mnt/storage/data/local/mol_med/${data_type}/${exp_data_dir}/samples"

if [ $jointype == "raw" ]
then
  hit="Hit"
else
  hit="Gene"
fi

# FIXME more checks/output here
# TODO add || exit 1 with a proper error correction after the main commands
mkdir -p $outdir && cd $outdir
mkdir -p stats collisions full
for i in $indices
do
    if [[ $tool == "migec" ]]
    then
      cat ${data_prefix}/${resdir}/${R}${m}${bc}${ql}/cdrfinal_${i}/S0_${R}.csv | awk -F"\t" '{print $1"\t"$3"\t"$4"\t"$5"\t"$6"\t"$2}' > mid${i}_clones.csv
      cat ${data_prefix}/${resdir}/${R}${m}${bc}${ql}/cdrfinal_${i}/S0_${R}.csv > full/MID${i}.csv
    else
      mixcr exportClones -c $chain -t -o -count -nFeature CDR3 -aaFeature CDR3 -vHit -jHit -dHit -cHit -fraction -f ${data_prefix}/${resdir}/${R}${m}${bc}${ql}/mixcr_${i}/analysis.*clns mid${i}_clones.csv
      mixcr exportClones -c $chain -t -o -vGene -jGene -dGene -cGene -p full -f ${data_prefix}/${resdir}/${R}${m}${bc}${ql}/mixcr_${i}/analysis.*clns full/MID${i}.csv
    fi
    cat mid${i}_clones.csv |  awk '{print $2}' | sort -k1,1 | uniq -c | awk '{if($1 > 1){print $0}}' | sort -k1,1 -r > collisions/mid${i}_collision.csv
done
cp ${data_prefix}/${resdir}/${R}${m}${bc}${ql}/*.csv stats/
cp ${data_prefix}/r{1,2}_stats.csv stats/
if [[ -f "${data_prefix}/../mid_labels.csv" ]]
then
  cp ${data_prefix}/../mid_labels.csv .
elif [[ -f "${data_prefix}/../${exp_dir}.csv" ]]
then
  cp ${data_prefix}/../${exp_dir}.csv mid_labels.csv
fi

# TODO: need a script to remove duplicates in place and recalculate counts and frequencies (for cristoph's scripts)
# might only be needed for t-cells as b-cells will be corrected via the below script

if [ $data_type == "bcr" ]
then
  if [ -n "$downsample" ]
  then
    downsample="-d ${downsample}"
  fi
  if [ -n "$groupby" ]
  then
    groupby="-g ${groupby}"
  fi
  Rscript /mnt/storage/data/code/clones2groups.R $downsample $groupby

  for i in $indices
  do
    cat mid${i}_clones.csv |  awk '{print $2}' | sort -k1,1 | uniq -c | awk '{if($1 > 1){print $0}}' | sort -k1,1 -r > collisions/mid${i}_collision.csv
  done
fi

/mnt/storage/data/code/join_mids.py --unique --type $jointype $indices

# TODO: either full or not, must be a choice
Rscript /mnt/storage/data/code/family_usage_full.R

# TODO: check if script exists? this should be automatised better, too - metadata should be pulled from somewhere
# and loaded into the R scripts
Rscript ../code/${execscript}.R

if [ $data_type == "bcr" ]
then
  if [ -d "fasta_groups" ]
  then
    cd fasta_groups
    for d in *
    do
      if [ -d $d ]
      then
        cd $d
        echo $d
        gctree.sh ${organism}
        cd ..
      fi
    done
  fi
fi
