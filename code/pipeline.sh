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
#germline="--impute-germline-on-export"
#contig="--contig-assembly"

# safe="_safe"
safe="_safe_long"

if [ $chain == "BCR" ] || [ $chain == "IGH" ] || [ $chain == "IGK" ] || [ $chain == "IGL" ]
then
  echo "###### B-cells ######"
  jprimers="c-primers"
  region="VDJRegion"
elif [ $chain == "TCR" ] || [ $chain == "TRA" ] || [ $chain == "TRB" ]
then
  echo "###### T-cells ######"
  jprimers="j-primers"
  region="CDR3"
fi

mkdir -p R12_m${m}${safe}

cd R12_m${m}${safe}

echo "##### R12_m${m}${safe} #####"

# summary
echo -e "MID\tInitial-reads\tInitial-UMIs\tOverseqFactor\tFinal-UMIs\tFinal-Clonotypes" > summary.csv
for i in $MIDS
do
  echo "#### MID${i} ####"

  nocache migec Checkout --skip-undef -p $p -oute ../../../../barcode${safe}.txt ../../MID${i}_*_R1_001.fastq ../../MID${i}_*_R2_001.fastq checkout_${i}
  nocache migec Histogram -p $p checkout_${i} histogram_${i}
  (cd histogram_${i} && Rscript ../../../../../../histogram.R)

  # if m non specified, extract collision and overseq thresholds from histogram output
  if [ $m == "b" ]
  then
    m2=`cat histogram_${i}/estimates.txt | tail -n1 | cut -f5`
  else
    m2=$m # one overall
  fi # TODO could be interesting in the future to do a mean/median approach, which however requires checkout & histogram to be run first for all samples, before assemble

  nocache migec Assemble $cf -p $p -m ${m2} checkout_${i}/S0_R1.fastq checkout_${i}/S0_R2.fastq assembly_${i}/

  mkdir -p mixcr_${i}

  nocache mixcr analyze amplicon -b imgt -t $p $contig $germline --report mixcr_${i}/report.log --verbose -f --species $organism --starting-material rna --5-end v-primers --3-end $jprimers --adapters adapters-present --receptor-type ${chain} --region-of-interest $region --only-productive --align "-OreadsLayout=Collinear" --assemble "-OseparateByC=false" --assemble "-OseparateByV=true" --assemble "-OseparateByJ=true" --assemble "-OqualityAggregationType=Average" --assemble "-OmaxBadPointsPercent=0" --assemble "-OcloneClusteringParameters=null" assembly_${i}/S0_R1.t${m2}${cf2}.fastq assembly_${i}/S0_R2.t${m2}${cf2}.fastq mixcr_${i}/analysis

  echo -e "$i\t`cat histogram_${i}/estimates.txt | awk 'NR>1{print $3"\t"$4"\t"$5}'`\t`cat mixcr_${i}/analysis.clonotypes.${chain}.txt | awk 'NR>1{sum+=$2} END {print sum"\t"NR}'`" >> summary.csv
done

cat summary.csv
