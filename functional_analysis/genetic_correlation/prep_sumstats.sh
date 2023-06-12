#! /bin/bash

# Make sure the downloaded Neale UKBB summary stats are in the correct format
# for LD score regression genetic correlation
#
# NOTE: Make sure to load the `ldsc` environment before running this script.

SUMSTATS=$1

# Add rsID to the summary stats so munge_sumstats.py works:
Rscript genetic_correlation/pre_munge.R ${SUMSTATS}

# Munge sumstats
munge_sumstats.py --sumstats ${SUMSTATS%%.tsv.gz}.clean.tsv.gz --snp SNP --N-col n_samples --a1 A1 --a2 A2 --p P --frq MAF --signed-sumstats beta,0  --ignore tstat --keep-maf --merge-alleles data/ldsc/ref_files/w_hm3.snplist --chunksize 100000 --out ${SUMSTATS%%.*}

# Remove the intermediate pre-munge file:
rm ${SUMSTATS%%.tsv.gz}.clean.tsv.gz
