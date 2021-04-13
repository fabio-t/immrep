#!/bin/bash

## removed chimeric and pseudogenes
migec="migec2"
data_type="tcr"
safe="_safe"
R="R12"

# indices=`echo {42..56}`
# exp_dir="170126_vancomycin_2"
# exp_data_dir="170126_vancomycin"

# indices=`echo {1..11} {13..16}`
# exp_dir="170707_diet_2"
# exp_data_dir="170707_diet"

# indices=`echo {1..11} {13..16}`
# exp_dir="170915_capture_2"
# exp_data_dir="170915_mixed"

# indices=`echo 17 {19..30}`
# exp_dir="170915_vancomycin"
# exp_data_dir="170915_mixed"

# indices=`echo 31 {34..37} 40 41`
# exp_dir="170915_braunschweig"
# exp_data_dir="170915_mixed"

# indices=`echo {19..39}`
# exp_dir="170927_capture_2"
# exp_data_dir="170927_capture"

# indices=`echo {40..64}`
# exp_dir="171024_capture_2"
# exp_data_dir="171024_capture"

indices=`echo {30..32} {34..42}`
exp_dir="210412_prevotella_rep"
exp_data_dir="210412_prevotella_rep"
R="R12_m1"
migec="umis"

## FIXME: originally there was this, but inside the script (for more recent
## experiments) we also add --type gene. To evaluate why not for Moritz
# ../../../../code/join_mids.py --unique $indices

# load and run script
source ../../../code/script.sh
