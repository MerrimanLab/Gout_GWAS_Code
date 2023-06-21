#! /bin/bash

# Quick script to make a list of variants in EUR Major GWAS summary stats

OUTPUT=$1

cut -f3 /Volumes/archive/merrimanlab/major_gwas_paper_archive/results_meta_analysis/EUR/full/EUR_meta_full1_clean_rsid.nfiltered.biallelic.txt > ${OUTPUT}/sumstat_variant_list.txt
cut -f1-2 /Volumes/archive/merrimanlab/major_gwas_paper_archive/results_meta_analysis/EUR/full/EUR_meta_full1_clean_rsid.nfiltered.biallelic.txt | awk '{print $1":"$2"_"}' > ${OUTPUT}/sumstat_variant_list.cpid.txt

# Also make a list with SNP IDs and allelic values IDs
Rscript src/cojo_script/make_avid.R /Volumes/archive/merrimanlab/major_gwas_paper_archive/results_meta_analysis/EUR/full/EUR_meta_full1_clean_rsid.nfiltered.biallelic.txt ${OUTPUT}/sumstats_avid.txt

cut -f2 ${OUTPUT}/sumstats_avid.txt > ${OUTPUT}/sumstats_avid_only.txt

