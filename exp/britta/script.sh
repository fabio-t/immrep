#!/bin/bash

data_type="bcr"
resdir="umis"
safe="_safe"
R="R12"
tool="migec"

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

# m="b"
# tool="mixcr"
# indices=`echo {1..11} 13`
# exp_dir="200911_replicates"
# exp_data_dir="200911_mouse_SI_technical_replicates"

# m="b"
# tool="mixcr"
# indices=`echo 43 {45..55} {57..69}`
# exp_dir="210329_diversity_si"
# exp_data_dir="210329_diversity_si"

m="b"
tool="mixcr"
indices=`echo 30 32 35 37 38 39 45 57 58 61 63 67 68 73`
exp_dir="210426_diversity_di"
exp_data_dir="210426_diversity_dual_indices"

# m="b"
# tool="mixcr"
# indices=`echo {1..15}`
# exp_dir="200709_comparison"
# exp_data_dir="200709_comparison_IgA_liver_gut_2"

# m="b"
# tool="mixcr"
# indices=`echo {16..24}`
# exp_dir="201215_sc_SPFvsMM12"
# exp_data_dir="201215_sc_SPFvsMM12"

# load and run script
source ../../code/script.sh
