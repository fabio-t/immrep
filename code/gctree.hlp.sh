#!/bin/bash

f=`basename $PWD`
name=`echo $f | cut -d_ -f2,3 | sed 's/IGH//g'`
deduplicate f.fasta --idmapfile idmap.csv --id_abundances --root ${name:0:10} --abundance_file abund.csv --frame 1 > f.phyi
mkconfig --quick f.phyi dnapars > dnapars.cfg
rm outtree outfile
phylip dnapars < dnapars.cfg 2>&1 | tee dnapars.log
gctree infer --verbose --root ${name:0:10} --frame 1 --idlabel outfile abund.csv
~/data_dir/code/ete.py
