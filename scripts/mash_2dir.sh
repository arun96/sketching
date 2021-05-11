#!/bin/bash
MSH=$1
# /Users/arun/Desktop/JHU/mash/mash-OSX64-v2.2/mash
DIR1=$2
DIR2=$3
TRGT=$4
for i in $DIR1/*; do
  iNM=$(basename "${i}" ".${i##*.}")
  iFN=${i##*/}
  #echo $iNM
  for j in $DIR2/*; do
    jNM=$(basename "${j}" ".${j##*.}")
    jFN=${j##*/}
    #echo $jNM
    $MSH dist $DIR1/${iFN} $DIR2/${jFN} -k 21 -p 8 >> $TRGT
  done
done
