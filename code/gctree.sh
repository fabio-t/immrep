rm -rf *.aln.* *.root.* gctree_*

shopt -s extglob

organism=${1:-mouse}
echo $organism

for f in IGH*@([0-9].fasta);
do
    ~/data_dir/code/gctree.single.sh $f $organism
done
