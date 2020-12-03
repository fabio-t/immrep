#!/bin/bash

migec="migec2"
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

indices=`echo {1..11} 13`
exp_dir="200911_replicates"
exp_data_dir="200911_mouse_SI_technical_replicates"

# load and run script
source ../../code/script.sh
