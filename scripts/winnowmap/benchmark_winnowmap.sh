#!/bin/bash
DB=$1
RDS=$2
OUTPUT=$3
KMERFILE=$4
for i in $RDS/*; do
    NM=${i##*/}
    echo "Mapping ${NM%%.*}"
    # minimap2 -a --secondary=no -t 12 $DB $RDS/${i} > $OUTPUT/${NM%%.*}.sam
    winnowmap -W $KMERFILE --secondary=no -I 16g  -t 50 $DB ${i} > $OUTPUT/${NM%%.*}.paf
done
echo "Mapping Done!"