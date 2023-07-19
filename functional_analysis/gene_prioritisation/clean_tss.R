# Script to clean up FANTOM5 TSS bed file
library(dplyr)
library(tidyr)

data = read.table('data/gene_prioritisation/TSS_human.bed', sep = '\t', header = F, skip = 1, stringsAsFactors = F)

# From the TSS perdiction README/PDF, peaks are classified into non-TSS,
# leniant TSS, and strict TSS (strict/leniant based on threshold).
#
# Remove the non-TSS, defined as grey (rgb = 211, 211, 211):
data = data %>% filter(V9 != '211,211,211')

# Filter out weird "p@chr*" type rows since you can't map it back to a gene:
data = data %>% select(V1:V4, V6) %>% filter(!grepl('p@chr', V4))
colnames(data) = c('chr', 'tss_start', 'tss_end', 'description', 'strand')

# Split out the gene names and scores
data = data %>% separate(description, remove = F, convert = T, sep = ',0', into = c('gene_list', 'score')) %>% separate_rows(gene_list, sep = ',') %>% separate(gene_list, into = c('p_num', 'gene'), sep = '@') %>% select(chr:tss_end, p_num:strand)

# Group by gene name and pull out the most supported TSS:
data = data %>% filter(p_num == 'p1')

# Clean up data
chrs = paste('chr', c(1:22, 'X'), sep = '')
data = data %>% filter(chr %in% chrs)
data$chr = gsub('chr', '', data$chr)
data$chr = as.numeric(gsub('X', '23', data$chr))

data = data %>% select(chr:tss_end, gene:strand)

# Use UCSC gene data to add on ENSG IDs, since some of the gencode gene symbols
# don't match up properly with the FANTOM5 gene symbol (e.g. AARS vs. AARS1)

ensembl = read.table('data/gene_prioritisation/UCSC_GRCh37_Genes_UniqueList2021.txt', sep = '\t', header = T, stringsAsFactors = F)
ensembl = ensembl %>% select(Gene, ensemblGeneID)

data = left_join(data, ensembl, by = c('gene' = 'Gene'))

write.table(data, 'data/gene_prioritisation/fantom5_tss.txt', sep = '\t', col.names = T, row.names = F, quote = F)

