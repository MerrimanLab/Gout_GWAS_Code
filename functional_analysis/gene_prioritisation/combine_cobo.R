library(dplyr)
library(purrr)
library(tidyr)

dat_list = c('data/gene_prioritisation/cobo_deglist.LPS.txt', 'data/gene_prioritisation/cobo_deglist.MSU.txt', 'data/gene_prioritisation/cobo_deglist.LPS_MSU.txt')

dat = map(dat_list, ~ read.table(.x, sep = '\t', header = T, stringsAsFactors = F) %>% mutate(contrast = gsub('[5h]', '', contrast), contrast = gsub('-', '_', contrast)))

# Remove any genes that have been converted into dates (e.g. Mar-01)
for (i in 1:length(dat)) {
	dat[[i]] = dat[[i]] %>% filter(!grepl('Mar|Sep', gene))
}

comb = full_join(dat[[1]], dat[[2]], by = 'gene', suffix = c('', '.MSU'))
comb = full_join(comb, dat[[3]], by = 'gene', suffix = c('.LPS', '.LPS_MSU')) %>% select(!contains('contrast'))
comb = comb %>% distinct

write.table(comb, 'data/gene_prioritisation/cobo_deglist.combined.txt', sep = '\t', col.names = T, row.names = F, quote = F)

