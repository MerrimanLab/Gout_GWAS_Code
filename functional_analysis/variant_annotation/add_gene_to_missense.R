library(dplyr)

missense = read.table('results/missense/candidate_variants.missense_info.txt', sep ='\t', header = T, stringsAsFactors = F)
missense = missense %>% filter(missense)

gene = read.table('data/missense/missense_gene_info.txt', sep ='\t', comment.char = '', header = T, stringsAsFactors = F)

res = left_join(missense, gene, by = c('SNP' = 'ID')) %>% select(SNP:ld_ancestry, INFO) %>% distinct %>% arrange(CHR, BP)
colnames(res)[ncol(res)] = 'gene'

write.table(res, 'results/missense/missense_gene.txt', sep = '\t', row.names = F, col.names = T, quote = F)

