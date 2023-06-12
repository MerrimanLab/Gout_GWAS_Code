#! /bin/bash

# Munge UKBB urate data
Rscript genetic_correlation/conv_lowp_rg.R data/urate/ukbb_full.rsid.txt

munge_sumstats.py --sumstats data/urate/ukbb_full.rsid.lowp.txt --snp SNP --N-col N --a1 minor --a2 major --p P --frq MAF --signed-sumstats effect,0 --keep-maf --merge-alleles data/ldsc/ref_files/w_hm3.snplist --out data/urate/ukbb_full.rsid --chunksize 1000000
