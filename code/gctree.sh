rm -rf *.aln.* *.root.* gctree_*

for f in *.fasta
do
    echo $f
    f1=`basename $f .fasta`.root.fasta
    echo ${f1}
    f2=`basename $f .fasta`.aln.fasta
    echo ${f2}
    cp ${f} ${f1}
    vname=`echo $f | cut -d_ -f1`
    vgerm=`sed -n "/${vname}\*01|/,/>/p" ~/data_dir/exp/britta/germline_fasta/mouse.v.ext.fasta | sed '1d; $d'`
    jname=`echo $f | cut -d_ -f2`
    jgerm=`sed -n "/${jname}\*01|/,/>/p" ~/data_dir/exp/britta/germline_fasta/mouse.j.ext.fasta | sed '1d; $d'`
    cl=`echo $f | cut -d_ -f3`
    div="3"
    n=`expr $cl / $div`
    nrep=`printf -- 'N%.0s' $(seq 1 ${n})`
    name=`echo ${vname}_${jname} | sed 's/IGH//g'`
    echo -e ">${name:0:10}\n${vgerm}" >> ${f1}
    echo -e "${nrep}\n${jgerm}" >> ${f1}
    muscle -in ${f1} -out ${f2}
done

for f in *.aln.fasta
do
    echo $f
    f2=${f%%.aln.fasta}
    echo $f2
    mkdir -p gctree_${f2}
    cd gctree_${f2}
    cat ../${f} | sed 's/>.*|abundance=\([0-9]\+\)/>\1/' > f.fasta
    cp ../${f} f.orig.fasta
    name=`echo $f | cut -d_ -f1,2 | sed 's/IGH//g'`
    deduplicate f.fasta --idmapfile idmap.csv --id_abundances --root ${name:0:10} --abundance_file abund.csv --frame 1 > f.phyi
    mkconfig --quick f.phyi dnaml > dnaml.cfg
    rm outtree outfile
    phylip dnaml < dnaml.cfg 2>&1 | tee dnaml.log
    gctree infer --verbose outfile abund.csv --root ${name:0:10} --frame 1 --idlabel
    ~/data_dir/code/ete.py
    cd ..
done
