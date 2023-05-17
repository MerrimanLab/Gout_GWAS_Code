#! /bin/bash

# Small script to process GUGC data set

SUMMARY=${1}.txt
FREQ=${1}.freq

# Uppercase alleles:
awk '{$6 = toupper($6); $7 = toupper($7); print $0}' ${SUMMARY} | sed -e 's/X/23/g' -e 's/Y/24/g' | tr ' ' '\t' > ${SUMMARY%.*}_chrpos.txt

# Add MAF:
paste ${SUMMARY%.*}_chrpos.txt <(cut -f4 ${FREQ}) > ${SUMMARY%.*}_chrpos_freq.txt

