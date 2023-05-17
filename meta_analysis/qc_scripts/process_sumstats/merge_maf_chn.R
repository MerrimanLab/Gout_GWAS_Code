# Small R script to add in MAF information to the summary stats.

library(dplyr)
library(readr)
library(data.table)

dat = read_delim('data/summary/EAS/gout_4653cases_4599controls.summary', delim = ' ')
indels = read_delim('data/summary/EAS/gout_4653cases_4599controls_indels.summary', delim = ' ')
maf = read_delim('data/summary/EAS/gout_4653cases_4599controls.maf', delim = ' ')

# Change column names, so there's less risk of puzzling major/minor:
colnames(dat) = c("CHR", "SNP", "POS", "major", "minor", "P", "effect", "SE")
colnames(indels) = c("CHR", "SNP", "POS", "major", "minor", "P", "effect", "SE")
colnames(maf) = c("CHR", "SNP", "minor", "major", "MAF_case", "MAF_control", "MAF_combined", "N")

# Replace "<" characters (for deletion variants) in dat to "<DEL>":
ind_minor = which(dat$minor == '<')
ind_major = which(dat$major == '<')
dat[ind_minor, 'minor'] = "<DEL>"
dat[ind_major, 'major'] = "<DEL>"

# Replace indel variants in dat with the variants from indels:
dat[which(dat$SNP %in% indels$SNP),] = indels

# Variants are in the same order for both dat and maf:
which(dat$SNP != maf$SNP) %>% length # 0

# Join A1 and A2:
maf$A1_A2 = paste(maf$minor, maf$major, sep = '_')
dat$A1_A2 = paste(dat$minor, dat$major, sep = '_')

# Check for the alleles that need to be flipped:
which(dat$A1_A2 != maf$A1_A2) %>% length

# Some of the variants need to be flipped:
ind = which(dat$A1_A2 != maf$A1_A2)

tmp = maf$minor[ind]
maf$minor[ind] = maf$major[ind]
maf$major[ind] = tmp
maf$A1_A2 = paste(maf$minor, maf$major, sep = '_')
maf$MAF_case[ind] = 1 - maf$MAF_case[ind]
maf$MAF_control[ind] = 1 - maf$MAF_control[ind]
maf$MAF_combined[ind] = 1 - maf$MAF_combined[ind]

# Check the order of the SNPs again:
which(dat$A1_A2 != maf$A1_A2) %>% length # 0

# Add MAF info to summary stats:
final = bind_cols(dat[,c(1:8)], maf[,c('MAF_case', 'MAF_control', 'MAF_combined', 'N')])

fwrite(final, 'data/summary/EAS/chinese_gout.tsv', quote = F, sep = '\t', row.names = F, col.names = T)
