# Quick script to make FATHMM and ABC input files
library(dplyr)

# Note that `dat` contains a bunch of variants that were flagged as missense,
# but haven't been checked whether it actually is a missense variant (which is
# done after candidate_variants.missense_info.txt is created)
dat = read.table('results/missense/candidate_variants.missense_info.txt', header = T, sep = '\t', stringsAsFactors = F)

# Load the final missense table so we know for sure which variants to remove as
# missense:
missense = read.table('results/missense/final_missense_table.txt', sep = '\t', header = T, stringsAsFactors = F)
miss_var = unique(missense$SNP)

# First remove all of the proxy variants to get back to the 1466 candidate variants:
dat = dat %>% filter(is.na(proxy_of))

# Then use the missense list to remove the 41 actual missense variants:
dat = dat %>% filter(!(SNP %in% miss_var))

# You can see that 4 variants flagged as missense, but isn't actually
# a missense when you look at the alleles in the GWAS
table(dat$missense)

# Create FATHMM input
fathmm_input1 = paste(dat$CHR, dat$BP, dat$minor, dat$major, sep = ',')
fathmm_input2 = paste(dat$CHR, dat$BP, dat$major, dat$minor, sep = ',')

fathmm_input = c(fathmm_input1, fathmm_input2)

writeLines(fathmm_input, 'results/missense/candidate_variants.fathmm_input.txt')

# Also make an input for ABC enhancer overlap (just changing the header)
abc = dat %>% select(SNP:BP)
colnames(abc)[3] = 'POS'
write.table(abc, 'results/missense/candidate_variants.abc_input.txt', sep = '\t', col.names = T, row.names = F, quote = F)

