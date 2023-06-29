#! /bin/bash

# Lazy script to edit the template file for DEPICT
# bash src/depict_scripts/edit_template.sh <data.txt> <output_prefix> <ancestry> <threshold> <value>

FILENAME=${PWD}/${1}
OUTPREFIX=$2
ANCESTRY=$3
THRES=$4
PCOL=$5

sed -e "s,DATA,$FILENAME,g" -e "s,OUTNAME,$OUTPREFIX,g" -e "s,THRESHOLD,$THRES,g" -e "s,VALUE,$PCOL,g" src/depict_scripts/template.${ANCESTRY}.cfg > src/depict_scripts/${OUTPREFIX}.cfg

