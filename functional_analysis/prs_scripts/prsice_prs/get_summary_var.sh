#! /bin/bash
#################################################################################

cut -f1 /Volumes/archive/merrimanlab/major_gwas_paper_archive/results_meta_analysis/EUR/full/EUR_meta_full.no_ukbb1.orig_snps.txt | tail -n+2 | tr '_' ' ' | awk '{printf("%02d:%d-%d\n", $1, $2, $2)}' > dat/prs/sumstats_varlist.txt
cut -f1 /Volumes/archive/merrimanlab/major_gwas_paper_archive/results_meta_analysis/EUR/male/EUR_meta_male.no_ukbb1.orig_snps.txt | tail -n+2 | tr '_' ' ' | awk '{printf("%02d:%d-%d\n", $1, $2, $2)}' >> dat/prs/sumstats_varlist.txt
cut -f1 /Volumes/archive/merrimanlab/major_gwas_paper_archive/results_meta_analysis/EUR/female/EUR_meta_female.no_ukbb1.orig_snps.txt | tail -n+2 | tr '_' ' ' | awk '{printf("%02d:%d-%d\n", $1, $2, $2)}' >> dat/prs/sumstats_varlist.txt

for i in {01..22}; do grep "^${i}" dat/prs/sumstats_varlist.txt | sort -V | uniq | tr '\n' ' ' > dat/prs/chr${i}_varlist.txt; done
