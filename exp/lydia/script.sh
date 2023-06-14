#!/bin/bash

data_type="tcr"
chain="TRA" # some experiments might require something different
resdir="umis"
bc="_safe_long"
R="R12"
tool="mixcr"
jointype="raw" # stripping them upstream now (8/12/2022)

# migec="migec"
# bc="_safe"

# indices=`echo {1..11} {13..17} {19..21}`
# exp_dir="170714_cultures"
# exp_data_dir="170714_cultures"

# indices=`echo {1..11} {13..17} {19..32}`
# exp_dir="171215_cultures"
# exp_data_dir="171215_cultures"

# indices=`echo {33..62}`
# exp_dir="171215_cultures_TRB"
# exp_data_dir="171215_cultures"

# indices=`echo {30..64}`
# exp_dir="181204_treg"
# exp_data_dir="181204_treg"

# indices=`echo 10 11 {13..17} 19 20 {22..24} {50..54} {62..65}`
# exp_dir="190218_irf4flox"
# exp_data_dir="190218_irf4flox"

# indices=`echo {24..46}`
# exp_dir="170912_capture"
# exp_data_dir="190320_capture"

# ## removed chimeric and pseudogenes
# indices=`echo {24..46}`
# exp_dir="170912_capture_2"
# exp_data_dir="190320_capture"
# migec="migec2"

# indices=`echo 48 50 51 52 54 55 56 58 59 57 60 62 63 61 49 53`
# exp_dir="190107_capture"
# exp_data_dir="190320_capture"

# indices=`echo 2 1 {4..11} 13 14 {25..33} {15..17} {19..24}`
# exp_dir="180830_treg"
# exp_data_dir="190715_treg"

# indices=`echo {35..44}`
# exp_dir="180830_sort"
# exp_data_dir="190715_treg"

# indices=`echo {20..55}`
# exp_dir="190220_treg"
# exp_data_dir="191104_mixed"

# indices=`echo {56..64}`
# exp_dir="190220_pools"
# exp_data_dir="191104_mixed"

# indices=`echo {27..38}`
# exp_dir="191213_itregtcon2"
# exp_data_dir="191213_mixed"

# indices=`echo {39..47}`
# exp_dir="190930_treg"
# exp_data_dir="191213_mixed"

# indices=`echo 55 {48..50} 52 53 54 51`
# exp_dir="191213_klrg1"
# exp_data_dir="191213_mixed"

# indices=`echo {20..55}`
# exp_dir="200415_resident"
# exp_data_dir="200415_resident"

# indices=`echo {1..11} {13..17} {19..47}`
# exp_dir="200615_diversity"
# exp_data_dir="200615_diversity"

## "strict" version, it requires a change in the $safe and $migec variables
# indices=`echo {35..44}`
# exp_dir="180830_sort_strict"
# exp_data_dir="190715_treg"
# migec="migec2"
# safe="_safe_long_2"

## removed chimeric and pseudogenes
# indices=`echo {1..11} {13..17} {19..47}`
# exp_dir="200615_diversity_2"
# exp_data_dir="200615_diversity"
# migec="migec2"

## removed chimeric and pseudogenes
# indices=`echo 230 240 250 260 {1..11} {13..17} {19..30}`
# exp_dir="200610_sc_seq"
# exp_data_dir="200610_sc_seq"
# migec="migec2"

## removed chimeric and pseudogenes
# indices=`echo {35..44}`
# exp_dir="180830_sort_2"
# exp_data_dir="190715_treg"
# migec="migec2"

## removed chimeric and pseudogenes
# indices=`echo 2 1 {4..11} 13 14 {25..33} {15..17} {19..24}`
# exp_dir="180830_treg_2"
# exp_data_dir="190715_treg"
# migec="migec2"

## removed chimeric and pseudogenes
# indices=`echo {30..64}`
# exp_dir="181204_treg_2"
# exp_data_dir="181204_treg"
# migec="migec2"

## removed chimeric and pseudogenes
# indices=`echo {20..55}`
# exp_dir="190220_treg_2"
# exp_data_dir="191104_mixed"
# migec="migec2"

## removed chimeric and pseudogenes
# indices=`echo {27..38}`
# exp_dir="191213_itregtcon2_2"
# exp_data_dir="191213_mixed"
# migec="migec2"

## removed chimeric and pseudogenes
# indices=`echo {39..47}`
# exp_dir="190930_treg_2"
# exp_data_dir="191213_mixed"
# migec="migec2"

## removed chimeric and pseudogenes
# indices=`echo {1..11} {13..17} {19..32}`
# exp_dir="171215_cultures_2"
# exp_data_dir="171215_cultures"
# migec="migec2"

## removed chimeric and pseudogenes
# indices=`echo {33..62}`
# exp_dir="171215_cultures_TRB_2"
# exp_data_dir="171215_cultures"
# migec="migec2"

## removed chimeric and pseudogenes
# indices=`echo {1..11} {13..17} {19..22}`
# exp_dir="200908_itregtc_2"
# exp_data_dir="200908_itregtc"
# migec="migec2"

### new pipeline moving forward

# indices=`echo {20..36}`
# exp_dir="211203_itregtc3"
# exp_data_dir="220330_itregtc3_klrg1"
# m="3"
# ql="20"

# indices=`echo {37..48}`
# exp_dir="211203_klrg1"
# exp_data_dir="220330_itregtc3_klrg1"
# m="3"
# ql="20"

# indices=`echo {1..11} 13`
# exp_dir="220425_klrg1"
# exp_data_dir="220425_klrg1_cd4cre"
# m="3"
# ql="30"

###### TODO: rerun with newer mixcr
# indices=`echo {45..55} {57..65}`
# exp_dir="220425_cd4cre"
# exp_data_dir="220425_klrg1_cd4cre"
# m="3"
# ql="30"

### TODO: rerun with newer mixcr
# indices=`echo {1..11} 13`
# exp_dir="211203_220425_klrg1_merged"
# exp_data_dir="220330_220425_klrg1_merged"
# m="3"
# ql="30"

# indices=`echo {1..11} {13..17} {19..31}`
# exp_dir="220218_itregtc4"
# exp_data_dir="220218_itregtc4"
# m="3"
# ql="30"
# bc="_safe_long2"

# indices=`echo 30 2 31 4 5 6 7 8 9`
# exp_dir="221007_hh-treg-exp1"
# exp_data_dir="221007_ana"
# m="3"
# ql="30"
# bc="_safe_long2"

## removed 17 and 19 as requested by Lydia; 02/06/2023
indices=`echo 10 11 13 14 15 16 20 21 22 23`
exp_dir="221007_anakadda2"
exp_data_dir="221007_ana"
m="3"
ql="30"
bc="_safe_long2"

# indices=`echo {24..29}`
# exp_dir="221007_clonekinetics"
# exp_data_dir="221007_ana"
# m="3"
# ql="30"
# bc="_safe_long2"

# indices=`echo {1..11} {13..16}`
# exp_dir="220928_cd4cre"
# exp_data_dir="220928_cd4cre"
# m="3"
# ql="30"
# bc="_safe_long2"

# indices=`echo {1..11} {13..17} {19..34}`
# exp_dir="221021_nktbonn"
# exp_data_dir="221021_nktbonn"
# m="3"
# ql="30"
# bc="_safe_long2"
# chain="TCR"

# indices=`echo {1..11} 13`
# exp_dir="221209_hh2"
# exp_data_dir="221209_hh2_mismatch12"
# m="3"
# ql="30"
# bc="_safe_long2"

# indices=`echo {14..17} {19..40}`
# exp_dir="221209_mm1-mm2"
# exp_data_dir="221209_hh2_mismatch12"
# m="3"
# ql="30"
# bc="_safe_long2"

# indices=`echo {1..11} {13..17} {19..26}`
# exp_dir="221118_itreg-tc5"
# exp_data_dir="221118_itreg-tc5"
# m="3"
# ql="30"
# bc="_safe_long2"

# indices=`echo {1..11} {13..17} 19 20`
# exp_dir="230104_treg-challenge"
# exp_data_dir="230317_treg-challenge_HH-exp3"
# m="3"
# ql="30"
# bc="_safe_long2"

# indices=`echo {21..32}`
# exp_dir="230104_HH_exp3"
# exp_data_dir="230317_treg-challenge_HH-exp3"
# m="3"
# ql="30"
# bc="_safe_long2"

# indices=`echo 1 3 4 6 7 9 10 13 14 16 17 20`
# exp_dir="230104_treg-challenge_merged"
# exp_data_dir="230317_treg-challenge_merged"
# m="3"
# ql="30"
# bc="_safe_long2"

# load and run script
source ../../code/script.sh
