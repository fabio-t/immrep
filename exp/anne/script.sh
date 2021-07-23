#!/bin/bash

data_type="tcr"
resdir="umis"
safe="_safe_long"
R="R12"
tool="mixcr"

# m="b"
# indices=`echo {25..39} {5..9} {45..49} 40 21 55 56 59`
# exp_dir="190218_treg"
# exp_data_dir="190214_irf4flow-resident"

m="2"
indices=`echo {20..27} 37 {28..36} {38..55}`
exp_dir="200415_resident"
exp_data_dir="200415_resident"

# load and run script
source ../../code/script.sh
