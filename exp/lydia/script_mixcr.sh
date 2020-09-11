#!/bin/bash

data_type="tcr"
safe="_safe"
R="R12"

# indices=`echo {35..44}`
# exp_dir="180830_sort_noUMIs"
# exp_data_dir="190715_treg"

indices=`echo {35..44}`
exp_dir="180830_sort_mixcrOnConsensus"
exp_data_dir="190715_treg"

# indices=`echo {1..11} {13..17} {19..47}`
# exp_dir="200615_diversity_mixcrOnConsensus"
# exp_data_dir="200615_diversity"

data_prefix="../../../local/mol_med/${data_type}/${exp_data_dir}/samples"

library="--library imgt"

align_opts="-t 4 ${library} -OreadsLayout=Collinear"
# -OreadsLayout=Collinear (might be needed - migec assembly are set to collinear R1/R2)

assemble_opts="-t 4 -OcloneClusteringParameters=null"
# -OassemblingFeatures=VDJRegion (with this, we might not need to further force clone uniqueness with join_mids)

mkdir -p $exp_dir && cd $exp_dir
rm -rf logs
mkdir -p full logs stats collisions
for i in $indices
do
  # MID_R1="${data_prefix}/*MID${i}_*_R1_001${trim}.fastq"
  MID_R1="${data_prefix}/migec/${R}/assembly_${i}${safe}/S0_R1.t1.fastq"

  # MID_R2="${data_prefix}/*MID${i}_*_R2_001${trim}.fastq"
  MID_R2="${data_prefix}/migec/${R}/assembly_${i}${safe}/S0_R2.t1.fastq"

  mixcr analyze amplicon \
        -s mouse \
        --starting-material rna \
        --5-end v-primers \
        --3-end j-primers \
        --adapters adapters-present \
        --receptor-type TCR \
        --region-of-interest VDJRegion \
        --only-productive \
        --align "${align_opts}" \
        --assemble "${assemble_opts}" \
        --export "-p full" \
        $MID_R1 $MID_R2 logs/mid${i} || exit 1
  mixcr exportClones \
        -count \
        -nFeature CDR3 \
        -aaFeature CDR3 \
        -vHit \
        -jHit \
        -fraction \
        -f logs/mid${i}.clna mid${i}_clones.csv || exit 1
  cp logs/mid${i}.clonotypes.TRA.txt full/mid${i}_clones.csv || exit 1
done

../../../code/join_mids.py --unique --type gene $indices

Rscript ../../../code/family_usage.R

Rscript ../code/${exp_dir}.R
