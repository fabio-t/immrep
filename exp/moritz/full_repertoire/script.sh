#!/bin/bash

migec="migec"
data_type="tcr"
safe="_safe"
R="R12"

# indices=`echo {42..56}`
# exp_dir="170126_vancomycin"
# exp_data_dir="170126_vancomycin"

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

## FIXME: originally there was this, but inside the script (for more recent
## experiments) we also add --type gene. To evaluate why not for Moritz
# ../../../../code/join_mids.py --unique $indices

# load and run script
source ../../../code/script.sh
