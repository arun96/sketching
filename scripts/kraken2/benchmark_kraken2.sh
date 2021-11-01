#!/bin/bash
DB=$1
RDS=$2
OUTPUT=$3
for i in $RDS/*; do
    NM=${i##*/}
    echo "Classifying ${NM%%.*}"
    kraken2 --db $DB --threads 12 --output $OUTPUT/${NM%%.*}.log --report $OUTPUT/report.txt ${i}
done
echo "Kraken2 is Done!"
