#! /bin/bash

# Script to clean the meta_analysis output from METAL

OUT=${1%.*}

# Make extra columns for CHR and POS; uppercase the alleles; remove markers
# with significant HetISq; retain only the columns required for TAMA (ID, chr,
# pos, MAF, beta, se, alleles, sample size).
awk ' NR == 1 {print "CHR", "POS", $0}; NR > 1 && $12 <= 95 { $2 = toupper($2); $3 = toupper($3); split($1, chrpos, "_"); print chrpos[1], chrpos[2], $0 }; ' ${1} | grep -v "NA" | tr ' ' '\t' | tee ${OUT}.clean.tmp | cut -f1-6,10-11,18 > ${OUT}.tmp

# Rearrange columns and add presence indicator:
awk ' NR == 1 {$10 = "presence"; print $3, $1, $2, $4, $5, $10, $9, $6, $7, $8}; NR > 1 {$10 = 1; print $3, $1, $2, $4, $5, $10, $9, $6, $7, $8}' ${OUT}.tmp > ${OUT}.tmp2

# Sort markers based on chromosome and position:
cat <(head -1 ${OUT}.tmp2) <(tail -n+2 ${OUT}.tmp2 | sort -n -k2,2 -k 3,3) > ${OUT}_pretama.tsv
cat <(head -1 ${OUT}.clean.tmp) <(tail -n+2 ${OUT}.clean.tmp | sort -n -k1,1 -k 2,2) > ${OUT}_clean.tsv

rm ${OUT}.clean.tmp* ${OUT}.tmp*
