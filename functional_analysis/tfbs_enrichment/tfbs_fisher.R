# Script to calculate the p-value of whether a TF is binding to the list of
# colocalised CpG sites by chance or not

library(dplyr)
library(purrr)

source('src/tfbs_enrichment/functions.R')

# Load in all the unique CpG sites present in the GoDMC data
cpgs = read.table('data/tfbs_data/godmc_cpg_sites.txt', header = T, sep = '\t', stringsAsFactors = F)
cpg_chrpos = read.table('data/meqtl_data/cpg450k_chrpos.txt', header = T, sep = '\t', stringsAsFactors = F)

sig_cpgs = read.table('results/meqtl_results/combined/meqtl_coloc.pp0.8_n100.all.txt', header = T, sep = '\t', stringsAsFactors = F)
sig_cpgs = sig_cpgs %>% mutate(CHR = cpg.CHR, BP = cpg.BP)

# Load RELI ChIP info
chip_info = read.table('data/tfbs_data/RELI_data/ChIPseq.index', header = T, sep = '\t', stringsAsFactors = F)
chip_info = chip_info %>% filter(Species == 'human') %>% arrange(label)

# Combine all TF data into one
all_tf = chip_info %>% pull(label) %>% paste('data/tfbs_data/RELI_data/ChIP-seq/', ., sep = '') %>% map(., ~ read.table(.x, header = F, sep = '\t', stringsAsFactors = F)) %>% map(., ~ .x %>% mutate(V1 = gsub('chr', '', V1, ignore.case = T)))
all_tf = simplify_tf_dat(all_tf, chip_info)

################################################################################
# Count how many cpg sites in total each TF binds to each meQTL sites
# NOTE: This will take long, so do it once and save the results
# count_overlap = map_dfc(all_tf, ~ check_overlap(cpgs, .x) %>% as.data.frame)

# tf_count = apply(count_overlap, 2, sum)
# names(tf_count) = names(all_tf)

# write.table(data.frame(tf = names(tf_count), count = tf_count), 'results/tfbs_results/tf_overlap_count.all_cpg.txt', row.names = F, col.names = T, sep = '\t', quote = F)

# Read in the TF counts for all the cpg sites present in the GoDMC data:
tf_dat = read.table('results/tfbs_results/tf_overlap_count.all_cpg.txt', header = T, sep = '\t', stringsAsFactors = F)
tf_count = as.vector(tf_dat$count)
names(tf_count) = tf_dat$tf

################################################################################
# Count how many cpg sites each TF binds to colocalised meQTL sites
coloc_cpgs = sig_cpgs %>% select(CHR, BP, mqtl_cpg) %>% distinct
coloc_count = map_dfc(all_tf, ~ check_overlap(coloc_cpgs, .x) %>% as.data.frame)

rownames(coloc_count) = coloc_cpgs$mqtl_cpg
colnames(coloc_count) = names(all_tf)

write.table(coloc_count, 'results/tfbs_results/tf_meqtl_matrix.txt', sep = '\t', col.names = T, row.names = T, quote = F)

# Count up how many CpG sites the TF bound to
bind_count = apply(coloc_count, 2, sum)
names(bind_count) = names(all_tf)

# Generate 2 by 2 contingency tables:
tables = map2(bind_count, tf_count, ~ make_table(.x, .y, set_n = nrow(sig_cpgs)))

# Run Fisher's test:
fisher_test = map(tables, ~ test_table(.x, test = 'fisher', return_p = F))
fisher_or = map_dbl(fisher_test, ~ .x$estimate)
fisher_l95 = map_dbl(fisher_test, ~ .x$conf.int[1])
fisher_u95 = map_dbl(fisher_test, ~ .x$conf.int[2])
fisher_p = map_dbl(tables, ~ test_table(.x, test = 'fisher'))

# List up all the significant TFs
fisher_sig = fisher_p[which(fisher_p < (0.05 / length(fisher_p)))]
fisher_sig = data.frame(transcription_factor = names(fisher_sig), fisher_p = fisher_sig, n_sig_cpg = bind_count[names(fisher_sig)])

# Add a column with a list of all the meQTL CpGs the TF bound to
cpg_list = map(coloc_count, ~ rownames(coloc_count)[.x])
sig_cpg_list = cpg_list[fisher_sig$transcription_factor]
fisher_sig$bound_cpgs = map_chr(sig_cpg_list, ~ paste(.x, collapse = ';'))

# Save output
write.table(fisher_sig, 'results/tfbs_results/sig_meqtl_coloc_tfs.txt', col.names = T, row.names = F, sep = '\t', quote = F)

