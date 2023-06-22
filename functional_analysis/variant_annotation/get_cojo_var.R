library(dplyr)
library(purrr)
library(tidyr)

cojo_list = list.files('data/cojo_dat', pattern = '.txt', full.names = T)
dat = map_dfr(cojo_list, ~ read.table(.x, sep = '\t', header = T, stringsAsFactors = F))

loci_list = list.files('data/loci_list', pattern = 'indep', full.names = T)
loci = map_dfr(loci_list, ~ read.table(.x, sep = '\t', header = T, stringsAsFactors = F) %>% select(SNP))
lead_snp = unique(loci$SNP)

# Make a list of all the variants found in the conditional analysis
cojo_var = unique(dat$SNP)

# Remove lead variants:
cojo_var = cojo_var[which(!(cojo_var %in% lead_snp))]

# Save the list
writeLines(cojo_var, 'data/missense/cojo_vars.txt')

