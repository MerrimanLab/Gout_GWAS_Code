#! /bin/bash

# Script to merge ChrX information to autosomes and clean up some of the
# summary stats (logging OR, calculating overall MAF, etc.).

# First merge the overall gout summary stats (no X) with the MAF info:
Rscript src/merge_maf_chn.R

# Convert the male-only and chrX data into a format so it can be merged
# with the other data:
parallel 'Rscript src/convert_chn.R {}' ::: $(ls data/summary/EAS/gout_male*.summary data/summary/EAS/gout_4653cases_4599controls_chrX.meta)

# # Merge X chromosome data with the rest of the summary stats for full sample:
# cat data/summary/EAS/chinese_gout.tsv <(awk ' NR > 1 {print $1, $2, $3, $5, $4, $10, $17, $9, $6, $7, $16, $15}' data/summary/EAS/gout_4653cases_4599controls_chrX.meta.converted | tr ' ' '\t') > data/summary/EAS/chinese_tmp && mv data/summary/EAS/chinese_tmp data/summary/EAS/chinese_gout.tsv

# # Male-only data:
# cat data/summary/EAS/gout_male_4210cases_4566controls.summary.converted <(tail -n+2 data/summary/EAS/gout_male_4210cases_4566controls_chrX.summary.converted) | awk ' NR > 1 {print $1, $2, $3, $5, $4, $10, $17, $9, $6, $7, $16, $15}' | tr ' ' '\t' > data/summary/EAS/chinese_gout_male.tsv
mv data/summary/EAS/gout_male_4210cases_4566controls.summary.converted data/summary/EAS/chinese_gout_male.tsv

