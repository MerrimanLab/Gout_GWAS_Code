#! /bin/bash

export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OMP_NUM_THREADS=1
export GOTO_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

SNP=$(echo ${1##*/} | sed 's/.txt//g')
SEX=$2

# Since the "clean.txt" file is arranged specifically for PAINTOR/FINEMAP
# program and doesn't have relevant columns from the ".txt" file, pull out the
# variants from the "clean.txt from the raw file (from the FINEMAP data):
cut -d ' ' -f1 ${1%%.*}.clean.txt | tail -n+2 > data/slalom_data/${SEX}/${SNP}.rsid.list

cat <(echo rsid chromosome position allele2 allele1 maf beta se z p) <(grep -Fwf data/slalom_data/${SEX}/${SNP}.rsid.list ${1}) | tr -s ' ' '\t' > data/slalom_data/${SEX}/${SNP}.clean.txt

# SLALOM assumes allele2 as OTH/effect/minor allele and allele1 as REF allele.
# Need to make sure the REF allele is actually REF allele (done by checking the
# Neale's UKBB variant list)

Rscript src/slalom/check_ref_allele.R data/slalom_data/${SEX}/${SNP}.clean.txt

