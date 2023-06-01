# Small R script to add in rsID information from SNP tracker to the summary
# statistics:

library(readr)
library(data.table)
library(dplyr)

args = commandArgs(trailingOnly = T)

dat = read_tsv(args[1])
beta = read_tsv(args[2])

# Remove the SNP column from beta and rename the study columns:
beta = beta[,2:ncol(beta)]

study1 = beta$study1[1]
study2 = beta$study2[1]
study3 = beta$study3[1]
study4 = beta$study4[1]

beta = beta %>% select(!contains("study"))

new_column = gsub('1', paste('_', study1, sep = ''), colnames(beta))
new_column = gsub('2', paste('_', study2, sep = ''), new_column)
new_column = gsub('3', paste('_', study3, sep = ''), new_column)
new_column = gsub('4', paste('_', study4, sep = ''), new_column)

colnames(beta) = new_column

rsid = read_delim(args[3], delim = ' ')
colnames(rsid) = c('cpid', 'SNP_original')

# Combine TAMA results, beta, and rsID info:
tmp = left_join(dat, beta, by = c('CHR', 'POS', 'cpid', 'N_STUD'))
tmp = left_join(tmp, rsid, by = 'cpid')

# Swap the input SNP ID with rsID, but keep the input SNP ID if there is no
# rsID for that variant:
tmp_snp = tmp$SNP
tmp$SNP_original = if_else(is.na(tmp$SNP_original), tmp$SNP, tmp$SNP_original)
tmp$SNP = tmp$SNP_original
tmp$SNP_original = tmp_snp

# Calculate the overal effect size using weighted average of the per-population
# effect size.

# Make two matrices - one for effect sizes and another for sample sizes:
effect = as.matrix(tmp[,c('mean_effect_AFR','mean_effect_EAS','mean_effect_EUR','mean_effect_LAT')])
variance = (as.matrix(tmp[,c('effect_sd_AFR','effect_sd_EAS','effect_sd_EUR','effect_sd_LAT')]))^2
weights = 1 / variance
sum_weights = apply(weights, 1, sum)
# sample_size = as.matrix(tmp[,c('N_AFR','N_EAS','N_EUR','N_LAT')])
# sum_n = apply(sample_size, 1, sum)
# weights = sample_size / sum_n

weighted_effect = apply(effect * weights, 1, sum) / sum_weights

# we_v = apply((weights * (effect - weighted_effect)^2), 1, sum)
# we_sd = sqrt(we_v)
# tmp$weighted_effect = weighted_effect
# tmp$weighted_effect_sd = we_sd

we_se = sqrt(1 / sum_weights)

tmp$weighted_effect = weighted_effect
tmp$weighted_effect_sd = we_se

fwrite(tmp, file = args[1], quote = F, sep = '\t', row.names = F, col.names = T)

