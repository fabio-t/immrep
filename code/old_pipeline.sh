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
echo -e "MID\tUMIs-in-table\tClones\tUMIs-in-Clones\tInitial-reads\tInitial-UMIs\tOverseqFactor" > R12_m${m}/summary.csv

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
        nocache migec CdrBlast --log-file R12_m${m}/cdrblast_${i}${safe}/log_raw.txt --log-overwrite --unmapped-fastq-file R12_m${m}/cdrblast_${i}${safe}/unmpapped_raw.fastq --cdr3-fastq-file R12_m${m}/cdrblast_${i}${safe}/cdr3_raw.fastq -q 15 -p $p --all-alleles $segments $chimeric -R $chain -S $organism R12_m${m}/checkout_${i}${safe}/S0_R1.fastq R12_m${m}/checkout_${i}${safe}/S0_R2.fastq R12_m${m}/cdrblast_${i}${safe}/S0_R12_raw.txt
        nocache migec CdrBlast --log-file R12_m${m}/cdrblast_${i}${safe}/log_asm.txt --log-overwrite --unmapped-fastq-file R12_m${m}/cdrblast_${i}${safe}/unmpapped_asm.fastq --cdr3-fastq-file R12_m${m}/cdrblast_${i}${safe}/cdr3_asm.fastq -a --cdr3-umi-table R12_m${m}/cdrblast_${i}${safe}/umitable.csv -q 15 -p $p --all-alleles $segments $chimeric -R $chain -S $organism R12_m${m}/assembly_${i}${safe}/S0_R1.t${m2}${cf2}.fastq R12_m${m}/assembly_${i}${safe}/S0_R2.t${m2}${cf2}.fastq R12_m${m}/cdrblast_${i}${safe}/S0_R12_asm.txt
        nocache migec FilterCdrBlastResults -p $p -n R12_m${m}/cdrblast_${i}${safe}/S0_R12_asm.txt R12_m${m}/cdrblast_${i}${safe}/S0_R12_raw.txt R12_m${m}/cdrfinal_${i}${safe}/S0_R12.csv
        nocache migec CreateCdrHypermGraph R12_m${m}/cdrfinal_${i}${safe}/S0_R12.csv R12_m${m}/net_${i}${safe}

        mkdir -p R12_m${m}/mixcr_${i}
        
        nocache mixcr analyze amplicon -t $p $contig $germline --report R12_m${m}/mixcr_${i}/report.log --verbose -f -b imgt --species $organism --starting-material rna --5-end v-primers --3-end c-primers --adapters adapters-present --receptor-type ${chain} --region-of-interest VDJRegion --only-productive --align "-OreadsLayout=Collinear" --assemble "-OseparateByC=false" --assemble "-OseparateByV=true" --assemble "-OseparateByJ=true" --assemble "-OqualityAggregationType=Max" --assemble "-OclusteringFilter.specificMutationProbability=1E-5" --assemble "-OmaxBadPointsPercent=0" R12_m${m}/assembly_${i}_safe/S0_R1.t${m2}${cf2}.fastq R12_m${m}/assembly_${i}_safe/S0_R2.t${m2}${cf2}.fastq R12_m${m}/mixcr_${i}/analysis

        echo -e "$i\t`tail -n+2 R12_m${m}/cdrblast_${i}${safe}/umitable.csv | wc -l`\t`tail -n+2 R12_m${m}/cdrfinal_${i}${safe}/S0_R12.csv | wc -l`\t`cat R12_m${m}/cdrfinal_${i}${safe}/S0_R12.csv | awk 'NR>1{s+=$1}END{print s}'`\t`cat R12_m${m}/histogram_${i}${safe}/estimates.txt | awk 'NR>1{print $3"\t"$4"\t"$5}'`" >> R12_m${m}/summary.csv

    done
done
