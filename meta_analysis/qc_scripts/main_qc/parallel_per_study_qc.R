################################################################################
# Parallel version of per_study_qc.R (use it with GNU parallel)
################################################################################

library(vroom)
library(dplyr)

options(scipen = 999)

args = commandArgs(trailingOnly = T)

# Load data:
dat = vroom(args[1])

dat_names = gsub('_clean.*', '', args[1])
dat_names = gsub('.*/', '', dat_names)

# Round sample size to closest integer, generate CHR:POS ID, add "allelic
# value" and generate MAC column:
dat$N = round(dat$N)
dat$MAC = 2 * pmin(dat$MAF, 1 - dat$MAF) * dat$N
dat$cpid = paste(dat$CHR, as.numeric(dat$BP), sep = '_')
minor_value = unlist(lapply(dat$minor, function(x) switch(x, A = 1, C = 2, G = 3, T = 5, -100)))
major_value = unlist(lapply(dat$major, function(x) switch(x, A = 1, C = 2, G = 3, T = 5, -100)))
dat$allelic_value = minor_value + major_value

################################################################################
# Per-study QC stuff
################################################################################

# Initialise a list to store some useful information:
bad_snp_list = c()
extreme_snp_list = c()
summary_list = c()

# A quick loop to check for bad SNPs:

# Check for bad SNPs:

# Count and report "weird" alleles (e.g. non-ACGT, multi-allelic):
bad_allele = which(dat$allelic_value < 0)

# Count and report P > 1, P < 0, NA, and infinite values:
bad_p = which(is.na(dat$P) | is.infinite(dat$P) | dat$P > 1 | dat$P < 0)

# Count and report NA, and infinite values:
bad_effect = which(is.na(dat$effect) | is.infinite(dat$effect))

# Count and report SE < 0, SE == NA, SE == inf:
bad_se = which(dat$SE <= 0 | is.na(dat$SE) | is.infinite(dat$SE))

# Count and report MAF < 0, > 1, monomorphic:
bad_maf = which(dat$MAF == 0 | dat$MAF == 1 | dat$MAF < 0 | dat$MAF > 1 | is.na(dat$MAF) | is.infinite(dat$MAF))

# Count and report variant with N < (0.05 * max(N)):
bad_n = which(dat$N < (0.05 * max(dat$N)) | is.na(dat$N))

# Grab all of the violating SNPs:
ind_bad = c(bad_allele, bad_p, bad_effect, bad_se, bad_maf, bad_n) %>% unique

# Report bad SNPs:
weird_snps = matrix(NA, 2, 8)
colnames(weird_snps) = c('Alleles', 'P', 'effect', 'SE', 'MAF', 'MAC', 'N', 'NumUniq')
weird_snps[1,] = c(length(bad_allele), length(bad_p), length(bad_effect), length(bad_se), length(bad_maf), NA, length(bad_n), length(ind_bad))

# Remove bad SNPs if present:
if (length(ind_bad) > 0) {
	bad_snps = dat[ind_bad,]
	dat = dat[-(ind_bad),]
} else {
	bad_snps = paste("There were no bad SNPs in data set", dat_names)
}

# Count and report extreme values (these will be noted down, but not removed
# from the analysis yet):
extreme_p = which(dat$P == 0)
extreme_effect = which(abs(dat$effect) > 10)
extreme_maf = which(dat$MAF <= 0.001 | dat$MAF >= 0.999)
extreme_mac = which(dat$MAC < (0.01 * dat$N))
ind_extreme = c(extreme_p, extreme_effect, extreme_maf, extreme_mac) %>% unique

# Report extreme SNPs:
weird_snps[2,] = c(NA, length(extreme_p), length(extreme_effect), NA, length(extreme_maf), length(extreme_mac), NA, length(ind_extreme))
rownames(weird_snps) = c('Bad', 'Extreme')

# Grab all of the extreme SNPs:
if (length(ind_extreme) > 0) {
	extreme_snps = dat[ind_extreme,]
} else {
	extreme_snps = paste("There were no extreme SNPs in data set", dat_names)
}

# Save some report data
sample_group = ifelse(grepl("full", args[1]), "full/", ifelse(grepl("female", args[1]), "female/", "male/"))

out_dir = gsub('[^/]*$', '', args[1])
out_dir = gsub('data/summary', 'results/pre_meta_qc', out_dir)
out_dir = paste(out_dir, sample_group, sep = '')

saveRDS(bad_snps, file = paste0(out_dir, dat_names, "_bad_snp_list.rds"))
saveRDS(extreme_snps, file = paste0(out_dir, dat_names, "_extreme_snp_list.rds"))
saveRDS(weird_snps, file = paste0(out_dir, dat_names, "_summary_report.rds"))

# Save the temporarily QCed data:
file_name = gsub('\\..*', '_qc.tsv', args[1])
vroom_write(dat, file_name)

