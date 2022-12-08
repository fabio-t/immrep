#!/bin/bash

data_type="bcr"
chain="IGH" # some experiments might require something different
resdir="umis"
bc="_safe_long"
R="R12"
tool="mixcr"
jointype="raw" # because it's clonal trees, we already strip alleles (FIXME not anymore, 15/04/2022 dev immunarch)

# exp_dir="181019_igg"
# exp_data_dir="181019_igg"
# indices=`echo 16 17 19 20 {22..28} 38 {61..69}`
# bc="_safe"
# m="b"
# ql=20

# m="b"
# indices=`echo {40..60}`
# bc="_safe"
# exp_dir="181105_iga"
# exp_data_dir="181105_iga"
# ql=15

# exp_dir="200623-5_iga-igg"
# exp_data_dir="200623_IgA-200625-IgG__IBD"
# indices=`echo {2..7} 11 {13..17} {19..21} {22..40}`
# bc="_safe"
# m="b"
# ql="mix"
# downsample="T"
# groupby="patient"

# exp_dir="200623_iga"
# exp_data_dir="200623_IgA-200625-IgG__IBD"
# indices=`echo {2..7} 11 {13..17} {19..21}`
# bc="_safe"
# m="b"
# ql="mix"
# downsample="T"
# groupby="patient"

# exp_dir="200911_replicates"
# exp_data_dir="200911_mouse_SI_technical_replicates"
# indices=`echo {1..11} 13`
# m="b"
# ql=20
# downsample=1000
# groupby="mouse"

# exp_dir="210329_diversity_si"
# exp_data_dir="210329_diversity_si"
# indices=`echo 43 {45..55} {57..69}`
# m="b"
# ql=20

# exp_dir="210426_diversity_di"
# exp_data_dir="210426_diversity_dual_indices"
# indices=`echo 30 32 35 37 38 39 45 57 58 61 63 67 68 73`
# m="b"
# ql=15

# exp_dir="200709_comparison"
# exp_data_dir="200709_comparison_IgA_liver_gut"
# indices=`echo {30..44}`
# m="b"
# ql=15
# downsample="T"
# groupby="mouse"

# execscript="200709_comparison"
# exp_dir="200709_comparison_2"
# exp_data_dir="200709_comparison_IgA_liver_gut_2"
# indices=`echo {1..15}`
# m="b"
# ql=15

# execscript="200709_comparison"
# exp_dir="200709_comparison_2short"
# exp_data_dir="200709_comparison_IgA_liver_gut_2short"
# indices=`echo {1..15}`
# m="b"
# ql=20

# exp_dir="201215_sc_SPFvsMM12"
# exp_data_dir="201215_sc_SPFvsMM12"
# indices=`echo {16..24}`
# m="b"
# ql=15

# exp_dir="220408_colonised"
# exp_data_dir="220408_colonised"
# indices=`echo {20..39}`
# m="b"
# ql="15"
# downsample="T"
# groupby="colonisation,name"

# exp_dir="220516_colonised"
# exp_data_dir="220516_colonised"
# indices=`echo {20..39}`
# m="4"
# ql="30"
# #downsample="1000"
# bc="_safe_long2"

exp_dir="220701_etc_biopsies"
exp_data_dir="220701_etc_biopsies"
indices=`echo {1..11} {13..17} {19..26} {28..34}`
m="3"
ql="30"
downsample="F"
bc="_safe_long2"


# load and run script
source ../../code/script.sh
