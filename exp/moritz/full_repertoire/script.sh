#!/bin/bash

indices=`echo {42..56}`
exp_dir="170126_vancomycin"
exp_data_dir="170126_vancomycin"

# indices=`echo {1..11} {13..16}`
# exp_dir="170707_diet"
# exp_data_dir="170707_diet"

# indices=`echo {1..11} {13..16}`
# exp_dir="170915_capture"
# exp_data_dir="170915_mixed"

# indices=`echo 17 {19..30}`
# exp_dir="170915_vancomycin"
# exp_data_dir="170915_mixed"

# indices=`echo 31 {34..37} 40 41`
# exp_dir="170915_braunschweig"
# exp_data_dir="170915_mixed"

# indices=`echo {19..39}`
# exp_dir="170927_capture"
# exp_data_dir="170927_capture"

# indices=`echo {40..64}`
# exp_dir="171024_capture"
# exp_data_dir="171024_capture"

mkdir -p $exp_dir && cd $exp_dir

R=R12
# R=R1

safe="_safe"
# safe=""

for i in $indices
do
   cat ../../../../local/mol_med/tcr_cd3/${exp_data_dir}/samples/migec/${R}/cdrfinal_${i}${safe}/S0_R12.csv | awk -F"\t" '{print $1"\t"$3"\t"$4"\t"$5"\t"$6"\t"$2}' > mid${i}_clones.csv
done

# ../../../../code/join_mids.py $indices
../../../../code/join_mids.py --unique $indices

Rscript ../../../../code/family_usage_full.R

Rscript ../code/script_${exp_dir}.R
