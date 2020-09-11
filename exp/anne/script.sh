#!/bin/bash

migec="migec"
data_type="tcr"
safe="_safe"
R="R12"

# indices=`echo {25..39} {5..9} {45..49} 40 21 55 56 59`
# exp_dir="190218_treg"
# exp_data_dir="190218_irf4flox"

indices=`echo {20..27} 37 {28..36} {38..55}`
exp_dir="200415_resident"
exp_data_dir="200415_resident"

# load and run script
source ../../code/script.sh
