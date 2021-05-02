#!/bin/bash

migec="migec"
data_type="bcr"
safe="_safe"
R="R12"

# indices=`echo 16 17 19 20 {22..28} 38 {61..69}`
# exp_dir="181019_igg"
# exp_data_dir="181019_igg"

# indices=`echo {40..60}`
# exp_dir="181105_iga"
# exp_data_dir="181105_iga"

# indices=`echo {2..7} 11 {13..17} {19..21}`
# exp_dir="200623_iga"
# exp_data_dir="200623_IBD_IgA"

# indices=`echo {22..40} 64`
# exp_dir="200625_igg"
# exp_data_dir="200625_IBD_IgG"

# indices=`echo {2..7} 11 {13..17} {19..21} {22..40} 64`
# exp_dir="200623-5_iga-igg"
# exp_data_dir="200623_IgA-200625-IgG__IBD"

# m="10"
# indices=`echo {1..11} 13`
# exp_dir="200911_replicates_migec-m${m}"
# exp_data_dir="200911_mouse_SI_technical_replicates"
# migec="umis"
# R="R12"

# m="1"
# indices=`echo 43 {45..55} {57..69}`
# exp_dir="210329_diversity_si_migec-m${m}"
# exp_data_dir="210329_diversity_si"
# migec="umis"
# R="R12"

m="1"
indices=`echo 30 32 35 37 38 39 45 57 58 61 63 67 68 73`
exp_dir="210426_diversity_di_migec-m${m}"
exp_data_dir="210426_diversity_dual_indices"
migec="umis"
R="R12"

# load and run script
source ../../code/script.sh
