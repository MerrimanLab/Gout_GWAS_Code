# Script to test the enrichment of 71 CHIP genes in the prioritisation table
library(dplyr)
library(tidyr)

gene = read.table('/Volumes/scratch/merrimanlab/major_gwas_paper_scratch/rikutakei/GWAS_Functional/data/gene_prioritisation/Gencode_GRCh37_Genes_UniqueList2021.txt', sep = '\t', header = T, stringsAsFactors = F)

protein = gene %>% filter(Coding == 'proteincoding') %>% distinct(ensemblGeneID, Gene) %>% mutate(ENSG_short = gsub('\\..*', '', ensemblGeneID))

dat = read.table('/Volumes/scratch/merrimanlab/major_gwas_paper_scratch/rikutakei/Gout_GWAS_Code/functional_analysis/gene_prioritisation/gene_prioritisation.norm_score.unique_overall.txt', sep = '\t', header = T, stringsAsFactors = F)
dat = dat %>% select(locus, gene, ENSG_short, contains('norm'))

# Reduce the list to protein coding genes:
dat = dat %>% filter(ENSG_short %in% protein$ENSG_short)

# How many CHIP genes are in the prioritisation list:
chip_gene = read.table('chip_genes.txt', sep = '\t', header = T, stringsAsFactors = F)

length(which(chip_gene$ENSG %in% dat$ENSG_short)) # 19

# Make contingency table and run chi-square test:
a = length(which(chip_gene$ENSG %in% dat$ENSG_short))
b = nrow(dat) - a
c = nrow(chip_gene) - a
d = nrow(protein) - (a + b + c)

mat = matrix(c(a, b, c, d), byrow = TRUE, 2, 2)

chisq.test(mat)

# Save a table of CHIP genes with the prioritisation scores added:
res = gene %>% mutate(ENSG_short = gsub('\\..*', '', ensemblGeneID)) %>% select(ENSG_short, Chrom:End) %>% left_join(chip_gene, ., by = c('ENSG' = 'ENSG_short'))
res = dat %>% select(ENSG_short, contains('norm')) %>% left_join(res, ., by = c('ENSG' = 'ENSG_short'))

write.table(res, 'chip_gene.prioritisation_score.txt', sep = '\t', col.names = T, row.names = F, quote = F)

