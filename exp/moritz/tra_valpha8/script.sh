#!/bin/bash

# trim="_trimmed"

# qual=30
# qual=20
# qual=22
# qual=15

# indices=`echo {1..24}`
# exp_dir="151102_ql${qual}${trim}"
# exp_data_dir="151102_neomycin"

## rerun with 3.0.13 and latest imgt
# indices=`echo {1..21}`
# exp_dir="160210_ampicillin"
# exp_data_dir="160210_ampicillin"

## rerun with 3.0.13 and latest imgt
# indices=`echo {1..39}`
# exp_dir="151021_neomycin"
# exp_data_dir="151021_neomycin"

# indices=`echo {1..39}`
# exp_dir="160203_ql${qual}${trim}"
# exp_data_dir="160203_ampicillin"

## rerun with 3.0.13 and latest imgt
indices=`echo {1..31}`
exp_dir="160817_prevotella"
exp_data_dir="160817_prevotella"

# indices=`echo {1..11} {13..17} {19..29}`
# exp_dir="161004_ql${qual}${trim}"
# exp_data_dir="161004_prevotella"

## rerun with 3.0.13 and latest imgt
# indices=`echo {1..21}`
# exp_dir="160516_diet"
# exp_data_dir="160516_diet"

## rerun with 3.0.13 and latest imgt
# indices=`echo {1..11} {13..17} {19..29} {32..41}`
# exp_dir="170126_vancomycin"
# exp_data_dir="170126_vancomycin"

# indices=`echo 17 {19..48} {50..55}`
# exp_dir="170707_ql${qual}${trim}"
# exp_data_dir="170707_diet"

# indices=`echo {1..11} {13..17} {19..44}`
# exp_dir="170901_vancomycin_ql${qual}${trim}"
# exp_data_dir="170901_vancomycin"

# indices=`echo {22..40}`
# exp_dir="170714_cultures"
# exp_data_dir="170714_cultures/samples"

library="--library imgt"
# library=""

mkdir -p $exp_dir && cd $exp_dir
echo `pwd`

# mkdir align clns logs full
mkdir -p logs full

rm logs/* full/*

for i in $indices
do
    # mixcr align -t 4 $library -r logs/align_mid${i}.txt -c TRA -s mouse -OreadsLayout=Opposite -OvParameters.geneFeatureToAlign=VTranscript -p rna-seq ../../../../local/mol_med/tcr_cd3/${exp_data_dir}/*MID${i}_*_R*001${trim}.fastq align/mid${i}_align.vdjca -f || exit 1
    # mixcr assemble -r logs/assemble_mid${i}.txt -OassemblingFeatures=CDR3 -OqualityAggregationType=Max -OcloneClusteringParameters.clusteringFilter.specificMutationProbability="1E-4" -ObadQualityThreshold=${qual} -f align/mid${i}_align.vdjca clns/mid${i}_clones.clns || exit 1
    # mixcr exportClones -t -o -count -nFeature CDR3 -aaFeature CDR3 -vHit -jHit -fraction -f clns/mid${i}_clones.clns mid${i}_clones.csv || exit 1
    # mixcr exportClones -t -o -p full -f clns/mid${i}_clones.clns full/mid${i}_clones.csv || exit 1

    mixcr analyze amplicon \
        -s mouse \
        --starting-material rna \
        --5-end v-primers --3-end j-primers \
        --adapters adapters-present \
        --receptor-type TRA \
        --only-productive \
        --threads 4 \
        --library imgt \
        --export "-p full" \
        ../../../../local/mol_med/tcr/${exp_data_dir}/samples/*MID${i}_*_R{1,2}_001${trim}.fastq logs/mid${i} || exit 1
    mixcr exportClones -count -nFeature CDR3 -aaFeature CDR3 -vHit -jHit -fraction -f logs/mid${i}.clns mid${i}_clones.csv || exit 1
    cp logs/mid${i}.clonotypes.TRA.txt full/mid${i}_clones.csv || exit 1
done

../../../../code/join_mids.py --unique --type gene $indices

mv mids_counts.csv mids_counts.all.csv
grep "TRAV12" mids_counts.all.csv > mids_counts.csv

# ../../../../code/make_stats.sh $indices > stats.csv
