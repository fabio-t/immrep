#!/bin/bash

set -e

chain=$1
shift

organism=$1
shift

MIDS="$@"

echo "MIDS: $MIDS"

echo "demultiplexing io.."

m=1
p=3
#chimeric="--allow-chimeric"
#segments="--all-segments"

# mkdir -p R1 R12
mkdir -p R12

# summary
echo -e "MID\tUMIs-in-table\tClones\tUMIs-in-Clones\tInitial-reads\tInitial-UMIs" > R12/summary.csv

for i in $MIDS
do
    for safe in "_safe"
    # for safe in "_safe_long_2"
    do
        # R1 only
        # migec Checkout -oute -m 15:0.2:0.2 ../../../barcode${safe}.txt ../MID${i}_*_R1_001.fastq . R1/checkout_${i}${safe}
        # migec Histogram R1/checkout_${i}${safe} R1/histogram_${i}${safe}
        # (cd R1/histogram_${i}${safe} && Rscript ../../../../../../../../code/histogram.R)
        # # FIXME: extract collision and overseq thresholds from histogram output?
        # migec Assemble --log-file R1/assembly_${i}${safe}/log.txt --log-overwrite -m $m --alignment-details R1/checkout_${i}${safe}/S0_R0.fastq . R1/assembly_${i}${safe}/
        # migec CdrBlast --log-file R1/cdrblast_${i}${safe}/log_raw.txt --log-overwrite --unmapped-fastq-file R1/cdrblast_${i}${safe}/unmpapped_raw.fastq --cdr3-fastq-file R1/cdrblast_${i}${safe}/cdr3_raw.fastq -q 15 -p $p --all-alleles $segments $chimeric -R $chain -S MusMusculus R1/checkout_${i}${safe}/S0_R0.fastq R1/cdrblast_${i}${safe}/S0_R0_raw.txt
        # migec CdrBlast --log-file R1/cdrblast_${i}${safe}/log_asm.txt --log-overwrite --unmapped-fastq-file R1/cdrblast_${i}${safe}/unmpapped_asm.fastq --cdr3-fastq-file R1/cdrblast_${i}${safe}/cdr3_asm.fastq -a --cdr3-umi-table R1/cdrblast_${i}${safe}/umitable.csv -q 15 -p $p --all-alleles $segments $chimeric -R $chain -S MusMusculus R1/assembly_${i}${safe}/S0_R0.t${m}.fastq R1/cdrblast_${i}${safe}/S0_R0_asm.txt
        # migec FilterCdrBlastResults -n R1/cdrblast_${i}${safe}/S0_R0_asm.txt R1/cdrblast_${i}${safe}/S0_R0_raw.txt R1/cdrfinal_${i}${safe}/S0_R0.csv

        # R1 and R2
        nocache migec Checkout -p $p -oute -m 15:0.2:0.2 ../../../barcode${safe}.txt ../MID${i}_*_R1_001.fastq ../MID${i}_*_R2_001.fastq R12/checkout_${i}${safe}
        nocache migec Histogram -p $p R12/checkout_${i}${safe} R12/histogram_${i}${safe}
        (cd R12/histogram_${i}${safe} && Rscript ../../../../../../../../code/histogram.R)
        # FIXME: extract collision and overseq thresholds from histogram output?
        nocache migec Assemble -p $p -m $m R12/checkout_${i}${safe}/S0_R1.fastq R12/checkout_${i}${safe}/S0_R2.fastq R12/assembly_${i}${safe}/
        nocache migec CdrBlast --log-file R12/cdrblast_${i}${safe}/log_raw.txt --log-overwrite --unmapped-fastq-file R12/cdrblast_${i}${safe}/unmpapped_raw.fastq --cdr3-fastq-file R12/cdrblast_${i}${safe}/cdr3_raw.fastq -q 15 -p $p --all-alleles $segments $chimeric -R $chain -S $organism R12/checkout_${i}${safe}/S0_R1.fastq R12/checkout_${i}${safe}/S0_R2.fastq R12/cdrblast_${i}${safe}/S0_R12_raw.txt
        nocache migec CdrBlast --log-file R12/cdrblast_${i}${safe}/log_asm.txt --log-overwrite --unmapped-fastq-file R12/cdrblast_${i}${safe}/unmpapped_asm.fastq --cdr3-fastq-file R12/cdrblast_${i}${safe}/cdr3_asm.fastq -a --cdr3-umi-table R12/cdrblast_${i}${safe}/umitable.csv -q 15 -p $p --all-alleles $segments $chimeric -R $chain -S $organism R12/assembly_${i}${safe}/S0_R1.t${m}.fastq R12/assembly_${i}${safe}/S0_R2.t${m}.fastq R12/cdrblast_${i}${safe}/S0_R12_asm.txt
        nocache migec FilterCdrBlastResults -p $p -n R12/cdrblast_${i}${safe}/S0_R12_asm.txt R12/cdrblast_${i}${safe}/S0_R12_raw.txt R12/cdrfinal_${i}${safe}/S0_R12.csv

        echo -e "$i\t`tail -n+2 R12/cdrblast_${i}${safe}/umitable.csv | wc -l`\t`tail -n+2 R12/cdrfinal_${i}${safe}/S0_R12.csv | wc -l`\t`cat R12/cdrfinal_${i}${safe}/S0_R12.csv | awk 'NR>1{s+=$1}END{print s}'`\t`cat R12/histogram_${i}${safe}/estimates.txt | awk 'NR>1{print $3"\t"$4}'`" >> R12/summary.csv

    done
done
