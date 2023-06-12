# Add trait info, pull out 934 primary phenotypes, and save significant results
library(dplyr)
library(tidyr)

args = commandArgs(trailingOnly = T)

dat = read.table(args[1], sep = '\t', header = T, stringsAsFactors = F)

trait_info = read.table('data/neale_ukbb/ukb31063_ldsc_sumstat_manifest.tsv', sep = '\t', header = T, stringsAsFactors = F)
trait_info = trait_info %>% select(phenotype:description, source:dilute, is_primary_gwas:ldsc_h2_significance) %>% filter(is_primary_gwas)

# Combine info
res = left_join(dat, trait_info, by = c('p2' = 'phenotype')) %>% arrange(desc(abs(rg)))
colnames(res) = c( "trait1", "phenotype_code", "rg", "se", "z", "p", "h2_obs", "h2_obs_se", "h2_int", "h2_int_se", "gcov_int", "gcov_int_se", "trait2", "source", "sex", "variable_type", "is_v2", "has_v2_replacement", "dilute", "is_primary_gwas", "ldsc_confidence", "ldsc_h2_significance")
res = res %>% select(trait1, trait2:variable_type, rg:gcov_int_se, is_primary_gwas:ldsc_h2_significance, phenotype_code) %>% filter(source != 'finngen') %>% filter(!is.na(ldsc_h2_significance)) %>% filter(ldsc_h2_significance != 'nonsig')

# Remove duplicates by choosing "better quality" trait (lower p-value if tied)
dup = res$trait2[which(duplicated(res$trait2))]

res$score1 = lapply(res$ldsc_confidence, function(x) switch(x, none = 0, low = 1, medium = 2, high = 3))
res$score2 = lapply(res$ldsc_h2_significance, function(x) switch(x, nonsig = 0, nominal = 1, z4 = 2, z7 = 3, 0))
res$score = as.numeric(res$score1) + as.numeric(res$score2)

remove = c()

for (i in 1:length(dup)) {
	ind = which(res$trait2 %in% dup[i])
	tmp = res[ind,]
	if (tmp$score[1] == tmp$score[2]) {
		low_sig = which(tmp$p != min(tmp$p))
		remove = c(remove, ind[low_sig])
	} else {
		low_sig = which(tmp$score == min(tmp$score))
		remove = c(remove, ind[low_sig])
	}
}

res = res[-remove,] %>% select(!contains('score'))

# Pull out significant genetic correlations
res_sig = res %>% filter(p <= 0.05/934)

# Save both files
filename = gsub('tmp', 'txt', args[1])
filename_sig = gsub('tmp', 'sig.txt', args[1])

write.table(res, filename, sep = '\t', col.names = T, row.names = F, quote = F)
write.table(res_sig, filename_sig, sep = '\t', col.names = T, row.names = F, quote = F)

