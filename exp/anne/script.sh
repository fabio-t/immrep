#!/bin/bash

data_type="tcr"
chain="TRA"
resdir="umis"
bc="_safe_long"
R="R12"
tool="mixcr"
jointype="gene" # because it's clonal trees, we already strip alleles

m="b"
indices=`echo {25..39} {5..9} {45..49} 40 21 55 56 59`
exp_dir="190218_treg"
exp_data_dir="190214_irf4flow-resident"
ql=20

# m="b"
# indices=`echo {20..27} 37 {28..36} {38..55}`
# exp_dir="200415_resident"
# exp_data_dir="200415_resident"
# ql=20

# load and run script
source ../../code/script.sh
