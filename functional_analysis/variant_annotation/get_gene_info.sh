#! /bin/bash
# Quick script to pull out gene info for the missense variants

cat <(zgrep '^#CHR' data/dbsnp_missense/missense.vcf.gz) <(zgrep -Fwf results/missense/missense_variants.txt data/dbsnp_missense/missense.vcf.gz) | sed 's/RS.*GENEINFO=//g' | sed 's/:.*//g' > data/missense/missense_gene_info.txt

