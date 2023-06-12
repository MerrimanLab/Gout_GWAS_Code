# Take a look at the difference/similarity of male-/female-specific gout GWAS
# genetic correlation results
library(dplyr)

male = read.table('results/genetic_correlation/gout_male/rg_result.sig.txt', sep = '\t', header = T, stringsAsFactors = F)
female = read.table('results/genetic_correlation/gout_female/rg_result.sig.txt', sep = '\t', header = T, stringsAsFactors = F)
full = read.table('results/genetic_correlation/gout/rg_result.sig.txt', sep = '\t', header = T, stringsAsFactors = F)

full = full %>% select(trait2, sex, rg, se, p, ldsc_confidence, ldsc_h2_significance)
female = female %>% select(trait2, sex, rg, se, p)
male = male %>% select(trait2, sex, rg, se, p)

res = full_join(full, female, by = c('trait2', 'sex'), suffix = c('', '.female'))
res = full_join(res, male, by = c('trait2', 'sex'), suffix = c('.full', '.male'))
res = res %>% select(trait2:sex, contains('rg'), contains('se'), contains('p'), contains('ldsc'))

# Test whether the difference in genetic correlation between men/women are
# significantly different

# The default value of "n_blocks" (number of block jacknife blocks) in ldsc.py
# is 200 - this is used to calculate the SE for rg, so use this (minus 1) as
# degrees of freedom
res = res %>% mutate(t = (rg.female - rg.male) / ((se.female ^ 2) + (se.male ^ 2)) ^ 0.5, p = 2 * pt(-abs(t), df = 199))

# Pull out traits with higher genetic correlation in female-specific GWAS than
# in male-specific GWAS
hi_female = res %>% filter(rg.female >= rg.male) %>% arrange(p)

# Save results
write.table(res, 'results/genetic_correlation/gout_sex_comparison.txt', sep = '\t', col.names = T, row.names = F, quote = F)
write.table(hi_female, 'results/genetic_correlation/gout_sex_comparison.hi_female.txt', sep = '\t', col.names = T, row.names = F, quote = F)
