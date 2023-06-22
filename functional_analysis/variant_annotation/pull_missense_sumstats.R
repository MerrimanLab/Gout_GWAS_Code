library(vroom)
library(purrr)
library(dplyr)

args = commandArgs(trailingOnly = T)

cohort = args[1]

file_list = c( '/Volumes/archive/merrimanlab/major_gwas_paper_archive/results_meta_analysis/EUR/ZZZ/EUR_meta_ZZZ1_clean_rsid.nfiltered.biallelic.txt', '/Volumes/archive/merrimanlab/major_gwas_paper_archive/results_meta_analysis/AFR/ZZZ/AFR_meta_ZZZ1_clean_rsid.nfiltered.biallelic.txt', '/Volumes/archive/merrimanlab/major_gwas_paper_archive/results_meta_analysis/EAS/ZZZ/EAS_meta_ZZZ1_clean_rsid.nfiltered.biallelic.txt', '/Volumes/archive/merrimanlab/major_gwas_paper_archive/results_meta_analysis/LAT/ZZZ/LAT_meta_ZZZ1_clean_rsid.nfiltered.biallelic.txt', '/Volumes/archive/merrimanlab/major_gwas_paper_archive/results_meta_analysis/TAMA/ZZZ/tama_ZZZ_clean.nfiltered.biallelic.txt')
file_list = gsub('ZZZ', cohort, file_list)

missense_list = readLines('results/missense/missense_variants.txt')

data = map(file_list, ~ vroom(.x) %>% filter(SNP %in% missense_list))

outfile = paste('results/missense/missense_variants.', c('EUR', 'AFR', 'EAS', 'LAT', 'TAMA'), '_info.', cohort, '.txt', sep = '')
map2(data, outfile, ~ write.table(.x, .y, sep = '\t', col.names = T, row.names = F, quote = F) )
