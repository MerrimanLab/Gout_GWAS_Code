#! /bin/bash

# Pull out all CpG sites used in (cis-)GoDMC data

cut -f1 /Volumes/archive/merrimanlab/reference_files/GoDMC_meQTL/cleaned/cis/GoDMC_cis_meQTL.full.tsv | tail -n+2 | sort --parallel=5 | uniq > data/tfbs_data/tmp

cat <(head -1 data/meqtl_data/cpg450k_chrpos.txt) <(grep -Fwf data/tfbs_data/tmp data/meqtl_data/cpg450k_chrpos.txt) > data/tfbs_data/godmc_cpg_sites.txt

rm data/tfbs_data/tmp
