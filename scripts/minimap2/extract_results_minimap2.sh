#!/bin/bash
DIR=$1
for i in $DIR/*; do
    NM=${i##*/}
    echo ">${NM%%.*}:"
    # samtools view ${i} | cut -f3 | sort | uniq -c
    cat ${i} | cut -f6 | sort | uniq -c
done
