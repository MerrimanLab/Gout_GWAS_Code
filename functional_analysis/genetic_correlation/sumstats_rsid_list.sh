#! /bin/bash

# Map the Neale UKBB summary stats variant ID to rsID

zcat data/neale_ukbb/full/age.gwas.imputed_v3.both_sexes.tsv.gz | cut -f1 > data/neale_ukbb/tmp

cat <(echo CHR POS A1 A2 ID) <(awk 'NR > 1 {split($1, info, ":"); print info[1], info[2], info[3], info[4], $1}' data/neale_ukbb/tmp) | tr ' ' '\t' > data/neale_ukbb/sumstats_variants.txt

# Run snptracker to map chr/pos to rsID
bash /Volumes/scratch/merrimanlab/CoreTools/CurrentVersions/snptracker/snptrack.sh data/neale_ukbb/sumstats_variants.txt

# Merge the results together in R
Rscript src/genetic_correlation/merge_rsid.R

