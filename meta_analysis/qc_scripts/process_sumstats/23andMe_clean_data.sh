#!/bin/bash

# Small script to clean up 23andMe data columns
#
# Two outputs:
# 1. *_maf.dat - contains 23andMe data with relevant columns from the raw data
#    and MAF calculated from dosage
# 2. *af_info.dat.gz - contains only the allele frequency info of *_maf.dat file

FILE=$1
OUT_NAME=${FILE%.*}

# Split the allele column into A and B alleles, remove SNPs that failed QC
# (i.e. "pass" column), and remove D/I SNPs:
paste <(cut -f5-8 $FILE) <(cut -f8 $FILE | sed -e 's/s/A/g' -e 's/\/.//g') <(cut -f8 $FILE | sed -e 's/s/B/g' -e 's/.\///g') <(cut -f11-23,25-26 $FILE) | awk ' $4 !~ "D" {if ($21 != "N") print}' | tr ' ' '\t' > ${OUT_NAME}_tmp

# Calculate allele count for A and B alleles based on the dosage
# information, and then calculate MAF. Also, rename some chromosomes:
awk ' BEGIN {CONVFMT = "%.8g"} { if (NR == 1) { $22 = "tot.A.0"; $23 = "tot.A.1"; $24 = "tot.B.0"; $25 = "tot.B.1"; $26 = "MAF.A.0"; $27 = "MAF.A.1"; $28 = "MAF.B.0"; $29 = "MAF.B.1"; $30 = "MAF.A.tot"; $31 = "MAF.B.tot"; $32 = "N"; print; next; } else if ($20 == "I") { $22 = (2 * $10) * (1 - $11); $23 = (2 * $12) * (1 - $13); $24 = (2 * $10) * $11; $25 = (2 * $12) * $13; $32 = $10 + $11; } else { $22 = $15 + (2 * $14); $23 = $18 + (2 * $17); $24 = $15 + (2 * $16); $25 = $18 + (2 * $19); $32 = $14 + $15 + $16 + $17 + $18 + $19; } if (($23 + $25) > 0 && ($22 + $24) > 0) { $26 = ($22 / ($22 + $24)); $27 = ($23 / ($23 + $25)); $28 = ($24 / ($22 + $24)); $29 = ($25 / ($23 + $25)); $30 = (($22 + $23) / ($22 + $23 + $24 + $25)); $31 = (($24 + $25) / ($22 + $23 + $24 + $25)); } else { $26 = "NA"; $27 = "NA"; $28 = "NA"; $29 = "NA"; $30 = "NA"; $31 = "NA"; } print; } ' ${OUT_NAME}_tmp | sed -e 's/chr//g' -e 's/^X/23/g' -e 's/^Y/24/g' | tr ' ' '\t' | tee ${OUT_NAME}_maf_tmp | cut -f1-3,5-6,10-19,22- > ${OUT_NAME}_af_info.dat

# Reduce the number of columns:
cut -f1-9,10,12,28-29,31-32 ${OUT_NAME}_maf_tmp > ${OUT_NAME}_maf.dat

# Remove temporary data and gzip *af_info.dat:
rm ${OUT_NAME}*tmp
gzip -f ${OUT_NAME}_af_info.dat

