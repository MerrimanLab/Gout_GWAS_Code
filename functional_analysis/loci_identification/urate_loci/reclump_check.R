# Check the (re)clump results
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

# Only interested in the merged version of the loci
loci_list = read.table('dat/urate_clumping/all_loci.merged_r2_0.2.txt', sep = '\t', header = T, stringsAsFactors = F)

lead_clump = lead_clump %>% filter(loci %in% loci_list$simplified_name)
lead = lead %>% filter(loci %in% loci_list$simplified_name)

# Find those that don't match up:
res = left_join(lead, lead_clump) %>% mutate(match = lead == clump_lead)

# Make sure these are different to the initial list of no matches
old_no_match = read.table('res/urate_clumping/lead_not_match_clump.txt', sep = '\t', header = T, stringsAsFactors = F)

not_match = res %>% filter(!match) %>% filter(!(loci %in% old_no_match$loci))

# None of the no matches are new, so all the merged loci have consistent clump
# lead variant
nrow(not_match)

# Save the list - this is the loci list for Tin+UKBB urate
# NOTE: Still need to consider independent signals from different clumps within
# locus, but I will use the lead variant at each locus for the purpose of coloc
write.table(res, 'res/urate_clumping/tin_ukbb.loci_list.txt', sep = '\t', col.names = T, row.names = F, quote = F)
