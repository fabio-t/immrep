#!/bin/bash

f=`basename $PWD`
name=`echo $f | cut -d_ -f2,3 | sed 's/IGH//g'`
deduplicate f.fasta --idmapfile idmap.csv --id_abundances --root ${name:0:10} --abundance_file abund.csv --frame 1 > f.phyi
mkconfig --quick f.phyi dnaml > dnaml.cfg
rm outtree outfile
phylip dnaml < dnaml.cfg 2>&1 | tee dnaml.log
gctree infer --verbose outfile abund.csv --root ${name:0:10} --frame 1 --idlabel
~/data_dir/code/ete.py
