# Script to load in all the PAINTOR annotation LDSC results and check which
# annotations should be used to run PAINTOR

library(dplyr)
library(tidyr)
library(purrr)

# Load annotation results
files = list.files('results/ldsc/paintor_ldsc/', pattern = 'results.txt')
files = paste('results/ldsc/paintor_ldsc/', files, sep = '')

annot_res = map_dfr(files, ~ read.table(.x, sep = '\t', header = T, stringsAsFactors = F)) %>% arrange(Coefficient_P_value)

# Pull out those that passed Bonferroni corrected threshold
thres = 0.05 / nrow(annot_res)

sig_annot = annot_res %>% filter(Coefficient_P_value <= thres) %>% mutate(Name = gsub('data/ldsc/1kgp_ref_paintor/tmp_annot/', '', Name)) %>% separate(Name, into = c('annot_origin', 'annotation'), sep = '/')

# Save table:
write.table(sig_annot, 'results/ldsc/paintor_ldsc/sig_annot.txt', sep = '\t', col.names = T, row.names = F, quote = F)
