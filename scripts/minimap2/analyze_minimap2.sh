#!/bin/bash
EXP=$1
GNM=$2
ALGN=$3
./extract_contigs.sh $GNM > ${EXP}_contigs
./extract_results_minimap2.sh $ALGN > ${EXP}_results
python analyze_results.py ${EXP}_contigs ${EXP}_results
