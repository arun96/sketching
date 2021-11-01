#!/bin/bash
DIR=$1
for i in $DIR/*; do
    NM=${i##*/}
    echo "Sequences in ${NM%%.*}:"
    grep '>' ${i}
done
