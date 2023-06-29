#! /bin/bash

# Quick script to convert summary stats into format that DEPICT can understand

awk 'NR == 1 {$2 = "POS"; print $0}; NR > 1 {$19 = gensub("_", ":", "g", $20); print $0}' /Volumes/archive/merrimanlab/major_gwas_paper_archive/results_meta_analysis/EUR/full/EUR_meta_full1_clean_rsid.nfiltered.biallelic.txt | tr ' ' '\t' > data/depict/eur_full.txt
awk 'NR == 1 {$2 = "POS"; print $0}; NR > 1 {$19 = gensub("_", ":", "g", $20); print $0}' /Volumes/archive/merrimanlab/major_gwas_paper_archive/results_meta_analysis/EAS/full/EAS_meta_full1_clean_rsid.nfiltered.biallelic.txt | tr ' ' '\t' > data/depict/eas_full.txt
awk 'NR == 1 {$2 = "POS"; print $0}; NR > 1 {$19 = gensub("_", ":", "g", $20); print $0}' /Volumes/archive/merrimanlab/major_gwas_paper_archive/results_meta_analysis/AFR/full/AFR_meta_full1_clean_rsid.nfiltered.biallelic.txt | tr ' ' '\t' > data/depict/afr_full.txt
awk 'NR == 1 {$2 = "POS"; print $0}; NR > 1 {$19 = gensub("_", ":", "g", $20); print $0}' /Volumes/archive/merrimanlab/major_gwas_paper_archive/results_meta_analysis/LAT/full/LAT_meta_full1_clean_rsid.nfiltered.biallelic.txt | tr ' ' '\t' > data/depict/lat_full.txt

