# Script to check which cpg site passes Bonferroni-corrected p-value and
# whether they're in cis with the lead variant

library(dplyr)
library(purrr)
library(tidyr)

args = commandArgs(trailingOnly = T)

# Load CpG site chr/bp info
cpg_info = read.table('data/meqtl_data/cpg450k_chrpos.txt', sep = '\t', header = T, stringsAsFactors = F)

# List out all the files and get variant cpid
files = list.files(args[1], pattern = 'mQTL.txt', full.names = T)
variants = gsub('.*/', '', files)
variants = gsub('_mQTL.txt', '', variants)

# Load in all the files
full_list = map2_dfr(files, variants, ~ read.table(.x, sep = '\t', header = F, stringsAsFactors = F) %>% mutate(snp = .y))
colnames(full_list) = c('cpg', 'p', 'snp')

# Split out cpid into chr/bp
full_list = full_list %>% separate(snp, sep = '_', remove = F, into = c('CHR', 'BP'), convert = T)

# Filter for significant meQTLs:
threshold = 0.05 / nrow(full_list)

sig_list = full_list[which(full_list$p <= threshold), ]

# Add chr/bp info to the significant list
sig_list = left_join(sig_list, cpg_info, by = 'cpg', suffix = c('.gwas', '.cpg'))

# Make sure they are all in cis (<= 1Mb from GWAS lead)
sig_list = sig_list %>% mutate(BP_diff = abs(BP.gwas - BP.cpg)) %>% filter(BP_diff <= 1000000)

# Save the whole table as well as just the list of significant cpgs (for the
# next step in the pipeline)
full_table_name = paste(args[1], 'all_sig_mQTL.txt', sep = '/')
write.table(sig_list, full_table_name, sep = '\t', col.names = T, row.names = F, quote = F)

out_files = gsub('_mQTL', '_sig_mQTL', files)
out_list = map(variants, ~ sig_list %>% filter(snp == .x) %>% pull(cpg))

# Remove those with no significant cpgs after Bonferroni correction
ind = which(map(out_list, ~ length(.x)) == 0)

out_files = out_files[-ind]
out_list = out_list[-ind]

map2(out_files, out_list, ~ writeLines(.y, .x))

