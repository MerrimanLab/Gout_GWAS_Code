# Merge the full/male/female missense tables

library(dplyr)
library(tidyr)
library(purrr)

data_list = c( 'results/missense/missense_variants_information.full.txt', 'results/missense/missense_variants_information.male.txt', 'results/missense/missense_variants_information.female.txt')
data_list = map(data_list, ~ read.table(.x, sep = '\t', header = T, stringsAsFactors = F))

# Merge the locus info for the variants
cand = read.table('results/missense/candidate_variants.missense_info.txt', sep = '\t', header = T, stringsAsFactors = F)

data_list = map(data_list, ~ cand %>% distinct(locus, locus_ancestry, cohort, SNP, variant_from, proxy_of) %>% left_join(.x, ., by = c('SNP', 'proxy_of')) %>% select(locus:variant_from, SNP:CADD_PHRED))

# Filter for the variants from relevant cohorts:
data_list[[1]] = data_list[[1]] %>% filter(grepl('full', cohort)) %>% distinct
data_list[[2]] = data_list[[2]] %>% filter(grepl('\\<male\\>', cohort)) %>% filter(!grepl('full', cohort)) %>% distinct
data_list[[3]] = data_list[[3]] %>% filter(grepl('female', cohort)) %>% filter(!grepl('full|\\<male\\>', cohort)) %>% distinct

res = map_dfr(data_list, ~ .x)

# Number of missense variants that were not proxies:
res %>% filter(is.na(proxy_of)) %>% pull(SNP) %>% unique %>% length
res %>% filter(is.na(proxy_of)) %>% pull(SYMBOL) %>% unique %>% length

# Number of missense variants that were proxies of other variants:
res %>% filter(!is.na(proxy_of)) %>% pull(SNP) %>% unique %>% length
res %>% filter(!is.na(proxy_of)) %>% pull(SYMBOL) %>% unique %>% length

# Save output
write.table(res, 'results/missense/final_missense_table.txt', sep = '\t', col.names = T, row.names = F, quote = F)
