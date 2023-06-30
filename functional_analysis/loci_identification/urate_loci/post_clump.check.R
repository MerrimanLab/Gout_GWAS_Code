# Check the clump results and make sure the lead variant in the input data is
# also the lead variant of the first clump. (Needs to be done because there
# were some warnings from plink saying that lead variant wasn't present in the
# 1KGP data).
library(dplyr)
library(purrr)
library(tidyr)

# Load all clump files
file_list = list.files('res/urate_clumping/', pattern = 'chr.*.clumped', full.names = T)
clump_list = map(file_list, ~ read.delim(.x, sep = '', header = T, stringsAsFactors = F) %>% select(CHR, BP, SNP, P))

loci_names = gsub('.clumped', '', file_list)
loci_names = gsub('.*//', '', loci_names)

names(clump_list) = loci_names

# Get lead variant from each input locus:
lead_list = map(loci_names, ~ read.table(paste('dat/urate_clumping/', .x, '.txt', sep = ''), sep = '\t', header = T, stringsAsFactors = F) %>% filter(P == min(P)) %>% pull(SNP) %>% paste(., collapse = ';'))

lead = data.frame(loci = loci_names, lead = unlist(lead_list)) %>% separate_rows(lead, sep = ';')

# Get lead variant from clump result:
lead_clump = map(clump_list, ~ .x %>% pull(SNP) %>% .[1])
lead_clump = data.frame(loci = names(clump_list), clump_lead = unlist(lead_clump))

# Find those that don't match up:
res = left_join(lead, lead_clump) %>% mutate(match = lead == clump_lead)

not_match = res %>% filter(!match)

# For these variants, generate locuszoom to check if the lead variant is just
# an anomaly or not
write.table(not_match, 'res/urate_clumping/lead_not_match_clump.txt', sep = '\t', col.names = T, row.names = F, quote = F)

# Save the whole list for checking LD between the lead variants across locus
write.table(res, 'res/urate_clumping/all_clump_lead.txt', sep = '\t', col.names = T, row.names = F, quote = F)
