#! /bin/bash

# Small script to process Partners data set

SUMMARY=$1

# Split the SNP column to get CHR/POS information:
awk 'NR == 1 {print "CHR", "POS", $0}; NR > 1 {split($1, chrpos, "_"); print chrpos[1], chrpos[2], $0}' ${SUMMARY} | tr ' ' '\t' > ${SUMMARY%.*}_chrpos.txt

# Add MAF information:
Rscript src/add_maf_partners.R ${SUMMARY%.*}_chrpos.txt

