rm -rf *.aln.* *.root.* gctree_*

shopt -s extglob

organism=${1:-mouse}
echo $organism

for f in IGH*@([0-9].fasta);
do
    echo $f
    d=`basename $f .fasta`
    echo $d
    f1=${d}.root.fasta
    echo ${f1}
    cp ${f} ${f1}

    vname=`echo $f | cut -d_ -f1`
    vgerm=`sed -n "/${vname}\*01|/,/>/p" ~/data_dir/exp/britta/germline_fasta/${organism}.v.ext.fasta | sed '1d; $d'`
    jname=`echo $f | cut -d_ -f2`
    jgerm=`sed -n "/${jname}\*01|/,/>/p" ~/data_dir/exp/britta/germline_fasta/${organism}.j.ext.fasta | sed '1d; $d'`
    cl=`echo $f | cut -d_ -f3`
    div="1" # long N sequence seems to help with alignment
    n=`expr $cl / $div`
    nrep=`printf -- 'N%.0s' $(seq 1 ${n})`
    name=`echo ${vname}_${jname} | sed 's/IGH//g'`
    echo -e ">${name:0:10}\n${vgerm}" >> ${f1}
    echo -e "${nrep}" >> ${f1}
    echo -e "${jgerm}" >> ${f1}

    f2=${d}.aln.fasta
    echo ${f2}
    # muscle -in ${f1} -out ${f2}
    Rscript ~/data_dir/code/align.R ${f1} ${f2}

    mkdir -p gctree_${d}
    cd gctree_${d}
    cp ../${f2} f.orig.fasta
    cat f.orig.fasta | sed 's/>.*|abundance=\([0-9]\+\)/>\1/' > f.fasta
    name=`echo $f | cut -d_ -f1,2 | sed 's/IGH//g'`
    deduplicate f.fasta --idmapfile idmap.csv --id_abundances --root ${name:0:10} --abundance_file abund.csv --frame 1 > f.phyi
    mkconfig --quick f.phyi dnaml > dnaml.cfg
    rm outtree outfile
    phylip dnaml < dnaml.cfg 2>&1 | tee dnaml.log
    gctree infer --verbose outfile abund.csv --root ${name:0:10} --frame 1 --idlabel
    ~/data_dir/code/ete.py
    cd ..
done
