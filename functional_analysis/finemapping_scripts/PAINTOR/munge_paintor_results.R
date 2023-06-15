# Script to pull out relevant information from the SuSiE results
library(dplyr)
library(purrr)
library(stringr)
library(tidyr)
library(vroom)

args = commandArgs(trailingOnly = T)

# First argument should be the directory where the SuSiE outputs are located
dir = args[1]
sex = args[2]

# Generate a list of files in the output directory
files = list.files(dir, pattern = 'results', full.names = T)
files = files[2:length(files)]
snp = gsub('.*/', '', files)
snp = gsub('.results', '', snp)

# Load EUR loci list
loci_file = gsub('sex', sex, 'data/loci_list/sex_indepSNP_summary_1jun2022.txt')
loci = read.table(loci_file, sep = '\t', header = T, stringsAsFactors = F)
loci = loci %>% filter(EUR == 'EUR')

# Load EUR summary stats just for adding rsID later on
gwas_file = gsub('sex', sex, '/Volumes/archive/merrimanlab/major_gwas_paper_archive/results_meta_analysis/EUR/sex/EUR_meta_sex1_clean_rsid.nfiltered.biallelic.txt')
dat = vroom(gwas_file)
dat = dat %>% select(CHR:SNP, P)

# Load results from PAINTOR
res_list = map(files, ~ vroom(.x) %>% arrange(desc(Posterior_Prob)))
names(res_list) = snp

# Function to get 99% credible set by summing up variants until the total
# posterior probability reaches 99%
pull_cs = function(dat, cred_prob = 0.99) {
	tmp = dat %>% arrange(desc(Posterior_Prob))
	prob = 0
	ind = c()
	i = 1
	while (prob < cred_prob) {
		ind = c(ind, i)
		prob = prob + tmp$Posterior_Prob[i]
		i = i + 1
	}
	res = tmp[ind,]
	return(res)
}

# Set credible set probability cutoff.
# When PAINTOR is allowed to have >1 causal variant, the sum of posterior
# probability at the locus becomes greater than 1
cs_prob = map_dbl(res_list, ~ sum(.x$Posterior_Prob) * 0.99)

# Generate 99% credible set
cs_list = imap(res_list, ~ pull_cs(.x, cred_prob = cs_prob[.y]) %>% mutate(locus.SNP = .y))

# Count the number of variants in the 99% credible set at each locus
cred_num = map_dbl(cs_list, ~ nrow(.x)) %>% as.data.frame
cred_num$SNP = rownames(cred_num)
colnames(cred_num) = c('cred_num', 'SNP')
loci = left_join(loci, cred_num) %>% select(locus, CHR:SNP, SEX, EUR, cred_num)

# Add rsid and P info from summary stats:
full_cred = map_dfr(cs_list, ~ .x)
full_cred = left_join(full_cred, dat)
full_cred = full_cred %>% select(locus.SNP, SNP, CHR:Posterior_Prob, P)
colnames(full_cred)[11] = "P.gwas"

pip_0.1 = full_cred %>% filter(Posterior_Prob >= 0.1)
pip_0.5 = full_cred %>% filter(Posterior_Prob >= 0.5)

# Save into the same directory as where the input files are
out_name = gsub('sex', sex, 'credible_pip_list.sex.txt')
out_name = paste(dir, out_name, sep = '')

write.table(full_cred, out_name, sep = '\t', col.names = T, row.names = F, quote = F)
write.table(pip_0.1, gsub('txt', '0.1.txt', out_name), sep = '\t', col.names = T, row.names = F, quote = F)
write.table(pip_0.5, gsub('txt', '0.5.txt', out_name), sep = '\t', col.names = T, row.names = F, quote = F)
write.table(loci, gsub('_pip_list.*', paste(c('_summary', sex, 'txt'), collapse = '.'), out_name), sep = '\t', col.names = T, row.names = F, quote = F)

