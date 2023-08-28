#! /bin/bash

# Quick script to add gout phenotype to the UKBB fam file

# The fam file is sorted, so just paste the sorted phenotype file to it
paste <(cut -d ' ' -f1-5 dat/prs/ukb_geno/ukbb_chr1.fam) <(tail -n+2 src/prs/keep_ids.txt | sort -V | cut -f3) | tr '\t' ' ' > dat/prs/ukb_geno/ukbb_gout_pheno.fam

