#! /bin/bash

# Script to filter out variants with N < 0.2 * max(N) (or doesn't meet the
# N study requirement)

FILE=$1
OUT=${1%.*}

# Get SNP, P-value, and N from the summary stats:
cut -f3,12,13,18,22 ${FILE} | sed 's/_[ACGT]_[ACGT]//g' > ${OUT}.subset

# Remove/filter variants that have N < 0.2 * max(N), and make a list of
# significant variants:
Rscript filter_median.R ${OUT}.subset

# grep out these variants from the original file to make a filtered sum stats:
cut -f1 ${OUT}.nfilter > ${FILE%.*}.npass.list
grep -Fwf ${FILE%.*}.npass.list ${FILE} > ${FILE%.*}.nfiltered.txt

