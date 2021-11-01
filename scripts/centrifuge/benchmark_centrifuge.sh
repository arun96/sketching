#!/bin/bash
DB=$1
RDS=$2
OUTPUT=$3
for i in $RDS/*; do
    NM=${i##*/}
    echo "Classifying ${NM%%.*}"
    centrifuge -x $DB -f ${i} -S $OUTPUT/${NM%%.*}.log -p 12
    # centrifuge-kreport -x $DB -f ${i} -S $OUTPUT/${NM%%.*}.log -p 12
done
echo "Centrifuge is Done!"
