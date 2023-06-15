# Script to pull out relevant information from the FINEMAP results
library(dplyr)
library(purrr)
library(stringr)
library(tidyr)

args = commandArgs(trailingOnly = T)

# First argument should be the directory where the FINEMAP outputs are located
dir = args[1]
sex = args[2]
loci_list = args[3]

# Load EUR loci list
loci = read.table(loci_list, sep = '\t', header = T, stringsAsFactors = F)
loci = loci %>% filter(EUR == 'EUR')

# Load credible set info
snp = loci$SNP
file_name = paste(dir, snp, sep = '')
cred_name = paste(file_name, '.cred', sep = '')
config_name = paste(file_name, '.config', sep = '')

# Load causal configurations and pull out the number of variants in each top
# configuration:
config_list = map(config_name, ~ read.table(.x, sep = ' ', header = T, stringsAsFactors = F))
config_list = map(config_list, ~ .x %>% filter(config != 'null'))

# Generate a vector of correct .cred files:
cred_num = map_dbl(config_list, ~ .x$k[1])
cred_name = paste(cred_name, cred_num, sep = '')

cred_list = map(cred_name, ~ read.table(.x, sep = ' ', header = T, stringsAsFactors = F))

names(cred_list) = loci$SNP
names(config_list) = loci$SNP

# Function to change the cred_list into a long table format
make_long = function(x) {
	tmp = x %>% select(-index)
	res = tmp %>% pivot_longer(everything(),
							   names_to = c(".value", "set"),
							   names_pattern =   "(\\w+)([\\d+])",
							   values_drop_na = TRUE)
	res = res %>% rename(credible_set = set, SNP = cred, pip = prob) %>% arrange(credible_set)
	return(res)
}

# Change the list into long format, add the locus SNP to it, and rbind
# everything together
cred_list = map(cred_list, ~ make_long(.x))
cred_list = map2(cred_list, names(cred_list), ~ .x %>% mutate(total_cred_set = length(unique(.x$credible_set)), SNP.lead = .y))
cred_table = map_dfr(cred_list, ~ .x)

# Count number of variants in each credible set
cred_num = cred_table %>% group_by(SNP.lead, credible_set) %>% summarise(total_snp_in_cs = n())
cred_table = left_join(cred_table, cred_num)
cred_table = left_join(loci, cred_table, by = c('SNP' = 'SNP.lead'), suffix = c('.lead', '.credible')) %>% select(locus, total_cred_set, SNP, CHR:BP, SEX, EUR, credible_set, total_snp_in_cs, SNP.credible, pip)
colnames(cred_table)[3:5] = paste(colnames(cred_table)[3:5], 'lead', sep = '.')

# Summarise the information from the config file
config_info = loci %>% select(locus, CHR, BP, SNP)
config_info$causal_config_prob = map_dbl(config_list, ~ .x$prob[1])
config_info$causal_config = map_chr(config_list, ~ .x$config[1])

# Generate tables with all the variants with posterior inclusion
# probability (PIP) > 0.1 and > 0.5
pip_0.1 = cred_table %>% filter(pip >= 0.1)
pip_0.5 = cred_table %>% filter(pip >= 0.5)

# Save into the same directory as where the input files are
write.table(cred_table, paste(c(dir, 'credible_table.', sex, '.txt'), collapse = ''), sep = '\t', col.names = T, row.names = F, quote = F)
write.table(config_info, paste(c(dir, 'config_summary.', sex, '.txt'), collapse = ''), sep = '\t', col.names = T, row.names = F, quote = F)
write.table(pip_0.1, paste(c(dir, 'credible_pip_list.', sex, '.0.1.txt'), collapse = ''), sep = '\t', col.names = T, row.names = F, quote = F)
write.table(pip_0.5, paste(c(dir, 'credible_pip_list.', sex, '.0.5.txt'), collapse = ''), sep = '\t', col.names = T, row.names = F, quote = F)

