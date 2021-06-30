#!/bin/bash

set -e

# saner programming env: these switches turn some bugs into errors
#set -o errexit -o pipefail -o noclobber -o nounset

chain=$1
shift
echo "chain: $chain"

organism=$1
shift
echo "organism: $organism"

m=$1
shift
echo "m: $m"

MIDS="$@"
echo "MIDS: $MIDS"

echo "demultiplexing io.."

p=6
#chimeric="--allow-chimeric"
#segments="--all-segments"
cf="--filter-collisions"
cf2=".cf"
# germline="--impute-germline-on-export"
# contig="--contig-assembly"

mkdir -p R12_m${m}

# summary
echo -e "MID\tInitial-reads\tInitial-UMIs\tOverseqFactor\tFinal-UMIs\tFinal-Clonotypes" > R12_m${m}/summary.csv

for i in $MIDS
do
    for safe in "_safe"
    # for safe in "_safe_long_2"
    do
        echo "##### MID${i} ${safe} #####"

        nocache migec Checkout --skip-undef -p $p -oute ../../../barcode${safe}.txt ../MID${i}_*_R1_001.fastq ../MID${i}_*_R2_001.fastq R12_m${m}/checkout_${i}${safe}
        nocache migec Histogram -p $p R12_m${m}/checkout_${i}${safe} R12_m${m}/histogram_${i}${safe}
        (cd R12_m${m}/histogram_${i}${safe} && Rscript ../../../../../../../../code/histogram.R)

        # if m non specified, extract collision and overseq thresholds from histogram output
        if [ $m == "b" ]
        then
          m2=`cat R12_m${m}/histogram_${i}${safe}/estimates.txt | tail -n1 | cut -f5`
        else
          m2=$m # one overall
        fi # TODO could be interesting in the future to do a mean/median approach, which however requires checkout & historam to be run first, once

        nocache migec Assemble $cf -p $p -m ${m2} R12_m${m}/checkout_${i}${safe}/S0_R1.fastq R12_m${m}/checkout_${i}${safe}/S0_R2.fastq R12_m${m}/assembly_${i}${safe}/
        
        mkdir -p R12_m${m}/mixcr_${i}
        
        nocache mixcr analyze amplicon -t $p $contig $germline --report R12_m${m}/mixcr_${i}/report.log --verbose -f -b imgt --species $organism --starting-material rna --5-end v-primers --3-end c-primers --adapters adapters-present --receptor-type ${chain} --region-of-interest VDJRegion --only-productive --align "-OreadsLayout=Collinear" --assemble "-OseparateByC=false" --assemble "-OseparateByV=true" --assemble "-OseparateByJ=true" --assemble "-OqualityAggregationType=Average" --assemble "-OmaxBadPointsPercent=0" --assemble "-OcloneClusteringParameters=null"  R12_m${m}/assembly_${i}_safe/S0_R1.t${m2}${cf2}.fastq R12_m${m}/assembly_${i}_safe/S0_R2.t${m2}${cf2}.fastq R12_m${m}/mixcr_${i}/analysis
        
        echo -e "$i\t`cat R12_m${m}/histogram_${i}${safe}/estimates.txt | awk 'NR>1{print $3"\t"$4"\t"$5}'`\t`cat R12_m${m}/mixcr_${i}/analysis.clonotypes.${chain}.txt | awk 'NR>1{sum+=$2} END {print sum"\t"NR}'`" >> R12_m${m}/summary.csv

    done
done
