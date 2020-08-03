#!/bin/bash

set -e

usage()
{
	bs list project 77623556 -F Name -F DateCreated -F Project.Id -F TotalSize
	echo "Usage: $0 Project.Id EXP_DIR"
	exit 1
}

if [ $# -lt 2 ]
then
	usage
fi

pid=$1
exp_dir=$2

mkdir -p ${exp_dir}/samples && cd $exp_dir

for ds in `bs list datasets -F Name -F Project.Id -F Id --project-id=$pid --terse`
do
	echo $ds
	bs download datasets --id $ds --extension fastq.gz -o samples/${ds}
	tar xvfz samples/${ds}.tar.gz -C samples/
	rm samples/${ds}.tar.gz
done

rm samples/*.json; gunzip samples/*.gz
