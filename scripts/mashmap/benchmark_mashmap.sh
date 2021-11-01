#!/bin/bash
DB=$1
RDS=$2
OUTPUT=$3
for i in $RDS/*; do
    NM=${i##*/}
    echo "Mapping ${NM%%.*}"
    mashmap -t 16 -r $DB -q ${i} -o $OUTPUT/${NM%%.*}.log
done
echo "Mapping Done!"
