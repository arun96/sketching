#!/bin/bash
MSH=$1
# /Users/arun/Desktop/JHU/mash/mash-OSX64-v2.2/mash dist
DIR=$2
TRGT=$3
for i in $DIR/*; do
  iNM=$(basename "${i}" ".${i##*.}")
  iFN=${i##*/}
  #echo $iNM
  for j in $DIR/*; do
    jNM=$(basename "${j}" ".${j##*.}")
    jFN=${j##*/}
    #echo $jNM
    # $MSH $DIR/${iFN} $DIR/${jFN} -k 21 -p 8 > $TRGT/${iNM}_${jNM}.log
    $MSH dist $DIR/${iFN} $DIR/${jFN} -k 21 -p 8 >> $TRGT
  done
done
