rm -rf *.aln.* *.root.* gctree_IGHV*

for f in *.fasta
do
    echo $f
    f1=`basename $f .fasta`.root.fasta
    echo ${f1}
    f2=`basename $f .fasta`.aln.fasta
    echo ${f2}
    cp ${f} ${f1}
    vname=`echo $f | cut -d_ -f1`
    vgerm=`sed -n "/${vname}\*01|/,/>/p" ~/data_dir/exp/britta/germline_fasta/mouse.ext.fasta | sed '1d; $d'`
    echo -e ">${vname:0:10}\n${vgerm}" >> ${f1}
    muscle -in ${f1} -out ${f2}
done

for f in *.aln.fasta
do
    echo $f
    f2=${f%%.aln.fasta}
    echo $f2
    mkdir -p gctree_${f2}
    cd gctree_${f2}
    cat ../${f} | sed 's/>.*|\([0-9]\+\)/>\1/' > f.fasta
    vname=${f%%_*}
    deduplicate f.fasta --id_abundances --root ${vname:0:10} --abundance_file abund.csv --frame 1 > f.phyi
    mkconfig --quick f.phyi dnaml > dnaml.cfg
    rm outtree outfile
    phylip dnaml < dnaml.cfg 2>&1 | tee dnaml.log
    gctree infer --verbose outfile abund.csv --root ${vname:0:10} --frame 1
    cd ..
done

