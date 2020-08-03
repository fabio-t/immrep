#!/bin/bash

echo -n "TAGS"
for i in $@
do
    echo -n ",MID${i}"
done

echo

for tag in `cat tags.txt` NNNN
do
    echo -n "$tag"
    for i in $@
    do
        echo -n ","
        if [[ `grep $tag mid_$i/mid${i}_tags_final.txt -c` == 0 ]]; then
            echo -n "0"
        else
            grep $tag mid_${i}/mid${i}_tags_final.txt | awk '{print $1}' | tr -d "\n"
        fi
    done

    echo
done
