#!/bin/bash

shopt -s extglob

organism=${2:-mouse}
echo $organism

f=$1
echo $f
d=`basename $f .fasta`
echo $d

f1=${d}.root.fasta
echo ${f1}
cp ${f} ${f1}

vname=`echo $f | cut -d_ -f1`
vgerm=`cat ~/data_dir/exp/britta/germline_fasta/${organism}.v.ext.fasta | xargs -n2 | grep "${vname}\*01" | cut -d" " -f2`
jname=`echo $f | cut -d_ -f2`
jgerm=`cat ~/data_dir/exp/britta/germline_fasta/${organism}.v.ext.fasta | xargs -n2 | grep "${jname}\*01" | cut -d" " -f2`
cl=`echo $f | cut -d_ -f3`
div="1" # long N sequence seems to help with alignment
n=`expr $cl / $div`
nrep=`printf -- 'N%.0s' $(seq 1 ${n})`
name=`echo ${vname}_${jname} | sed 's/IGH//g'`
echo -e ">${name:0:10}\n${vgerm}" >> ${f1}
echo -e "${nrep}" >> ${f1}
echo -e "${jgerm}" >> ${f1}

rm -rf ${d}.aln.*
f2=${d}.aln.fasta
echo ${f2}
# muscle -in ${f1} -out ${f2}
Rscript ~/data_dir/code/align.R ${f1} ${f2}

rm -rf gctree_${d}
mkdir gctree_${d}
cd gctree_${d}
cp ../${f2} f.orig.fasta
cat f.orig.fasta | sed 's/>.*|abundance=\([0-9]\+\)/>\1/' > f.fasta
~/data_dir/code/gctree.hlp.sh
cd ..
