#!/bin/bash
JVD=$1
GNMS=$2
RDS=$3
SNP=$4
INS=$5
DEL=$6
MLN=$7
CVG=$8
OPT=$9
javac $JVD/*.java
for i in $GNMS/*; do
    if [ "${i}" != "${i%.fasta}" ] || [ "${i}" != "${i%.fna}" ] || [ "${i}" != "${i%.fa}" ];then
        GNM=${i##*/}
        echo "Generating Reads for ${GNM%%.*}"
        java -cp $JVD ReadSimulator ${i} $RDS/${GNM%%.*}.fasta $SNP $INS $DEL $MLN $CVG $OPT
    fi
done
echo "All Reads Generated!"
