# Script to convert TAMA logBF into something PLINK/DEPICT can understand
# (NOTE: This may not be necessary, if we decide to only look at EUR results)

# DEPICT uses PLINK to generate a list of loci using the summary stats P-value,
# so I'll have to convert logBF to some value that PLINK can understand.

library(vroom)
library(dplyr)

options(scipen = 99)

dat = vroom('/Volumes/archive/merrimanlab/major_gwas_paper_archive/results_meta_analysis/TAMA/full/tama_full_clean.nfiltered.biallelic.txt')

dat$value = ifelse(dat$logBF < 1, 1, 1 / dat$logBF)

dat = dat %>% select(SNP:OTH, value) %>% mutate(cpid = paste(CHR, POS, sep = ':'))

vroom_write(dat, 'data/depict/tama_full.txt')
