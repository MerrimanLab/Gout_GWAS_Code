# Merge rsID to the list of IDs in the sumstats
library(dplyr)

dat = read.table('data/neale_ukbb/sumstats_variants.txt', header = T, stringsAsFactors = F)
rsid = read.table('data/neale_ukbb/sumstats_variants.result.txt', header = T, stringsAsFactors = F)

dat$cpid = paste(dat$CHR, dat$POS, sep = '_')

res = left_join(dat, rsid, by = c('cpid' = 'CHR_BP')) %>% filter(!is.na(SNP)) %>% arrange(CHR, POS)

write.table(res, 'data/neale_ukbb/sumstats_variants.clean.txt', sep = '\t', col.names = T, row.names = F, quote = F)

