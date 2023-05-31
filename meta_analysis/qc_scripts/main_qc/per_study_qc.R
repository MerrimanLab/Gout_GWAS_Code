################################################################################
# Rscript for reporting the number of variants violating defined thresholds,
# and outputs a file containing the violating variants
################################################################################

library(data.table)
library(dplyr)

args = commandArgs(trailingOnly = T)

data_list = list()

dat_names = gsub('\\..*', '', args)
dat_names = gsub('.*/', '', dat_names)

# Load data:
for (i in 1:length(args)) {
	data_list[[i]] = as.data.frame(fread(args[i]))
}

names(data_list) = dat_names

# Round sample size to closest integer, generate CHR:POS ID, add "allelic
# value" and generate MAC column:
for (i in 1:length(data_list)) {
	tmp_dat = data_list[[i]]
	tmp_dat$N = round(tmp_dat$N)
	tmp_dat$MAC = 2 * pmin(tmp_dat$MAF, 1 - tmp_dat$MAF) * tmp_dat$N
	tmp_dat$cpid = paste(tmp_dat$CHR, tmp_dat$BP, sep = '_')
	minor_value = unlist(lapply(tmp_dat$minor, function(x) switch(x, A = 1, C = 2, G = 3, T = 5, -100)))
	major_value = unlist(lapply(tmp_dat$major, function(x) switch(x, A = 1, C = 2, G = 3, T = 5, -100)))
	tmp_dat$allelic_value = minor_value + major_value
	data_list[[i]] = tmp_dat
}

################################################################################
# Per-study QC stuff
################################################################################

# Initialise a list to store some useful information:
bad_snp_list = list()
extreme_snp_list = list()
summary_list = list()

# A quick loop to check for bad SNPs:
for (i in 1:length(data_list)) {
	tmp_dat = data_list[[i]]

	# Check for bad SNPs:

	# Count and report "weird" alleles (e.g. non-ACGT, multi-allelic):
	bad_allele = which(tmp_dat$allelic_value < 0)

	# Count and report P > 1, P < 0, NA, and infinite values:
	bad_p = which(is.na(tmp_dat$P) | is.infinite(tmp_dat$P) | tmp_dat$P > 1 | tmp_dat$P < 0)

	# Count and report NA, and infinite values:
	bad_effect = which(is.na(tmp_dat$effect) | is.infinite(tmp_dat$effect))

	# Count and report SE < 0, SE == NA, SE == inf:
	bad_se = which(tmp_dat$SE <= 0 | is.na(tmp_dat$SE) | is.infinite(tmp_dat$SE))

	# Count and report MAF < 0, > 1, monomorphic:
	bad_maf = which(tmp_dat$MAF == 0 | tmp_dat$MAF == 1 | tmp_dat$MAF < 0 | tmp_dat$MAF > 1 | is.na(tmp_dat$MAF) | is.infinite(tmp_dat$MAF))

	# Count and report variant with N < (0.05 * max(N)):
	bad_n = which(tmp_dat$N < (0.05 * max(tmp_dat$N)) | is.na(tmp_dat$N))

	# Grab all of the violating SNPs:
	ind_bad = c(bad_allele, bad_p, bad_effect, bad_se, bad_maf, bad_n) %>% unique

	# Report bad SNPs:
	weird_snps = matrix(NA, 2, 8)
	colnames(weird_snps) = c('Alleles', 'P', 'effect', 'SE', 'MAF', 'MAC', 'N', 'NumUniq')
	weird_snps[1,] = c(length(bad_allele), length(bad_p), length(bad_effect), length(bad_se), length(bad_maf), NA, length(bad_n), length(ind_bad))

	# Remove bad SNPs if present:
	if (length(ind_bad) > 0) {
		bad_snps = tmp_dat[ind_bad,]
		bad_snp_list[[i]] = bad_snps
		tmp_dat = tmp_dat[-(ind_bad),]
	} else {
		bad_snp_list[[i]] = paste("There were no bad SNPs in data set", i)
	}

	# Count and report extreme values (these will be noted down, but not removed
	# from the analysis):
	extreme_p = which(tmp_dat$P == 0 | tmp_dat$P < 1e-300)
	extreme_maf = which(tmp_dat$MAF <= 0.001 | tmp_dat$MAF >= 0.999)
	extreme_mac = which(tmp_dat$MAC < (0.01 * tmp_dat$N))
	extreme_n = which(tmp_dat$N < (0.01 * tmp_dat$N))
	ind_extreme = c(extreme_p, extreme_maf, extreme_mac, extreme_n) %>% unique

	# Report extreme SNPs:
	weird_snps[2,] = c(NA, length(extreme_p), NA, NA, length(extreme_maf), length(extreme_mac), length(extreme_n), length(ind_extreme))
	rownames(weird_snps) = c('Bad', 'Extreme')

	# Grab all of the extreme SNPs:
	if (length(ind_extreme) > 0) {
		extreme_snps = tmp_dat[ind_extreme,]
		extreme_snp_list[[i]] = extreme_snps
	} else {
		extreme_snp_list[[i]] = paste("There were no extreme SNPs in data set", i)
	}

	summary_list[[i]] = weird_snps
	data_list[[i]] = tmp_dat
}

# Now that the data is relatively clean, remove easy duplicates:

# Function to check duplicates and remove the counterpart with smaller MAC:
remove_dup = function(data) {
	dup = data$cpid[which(duplicated(data$cpid))] %>% unique
	if (length(dup) > 0) {
		snps_to_remove = c()
		dup_dat = data[data$cpid %in% dup, ]
		for(i in 1:length(dup)) {
			tmp = dup_dat[dup_dat$cpid %in% dup[i], ]
			snps_to_remove = c(snps_to_remove, rownames(tmp)[which(tmp$MAC != max(tmp$MAC))])
		}
		data = data[which(!(rownames(data) %in% snps_to_remove)), ]
	}
	return(data)
}

data_list = lapply(data_list, function(x) x = remove_dup(x))

# Save some report data
out_dir = gsub('[^/]*$', '', args[1])
out_dir = gsub('data/summary', 'results/pre_meta_qc', out_dir)

saveRDS(bad_snp_list, file = paste0(out_dir, "bad_snp_list.rds"))
saveRDS(extreme_snp_list, file = paste0(out_dir, "extreme_snp_list.rds"))
saveRDS(summary_list, file = paste0(out_dir, "summary_report.rds"))

# Save the temporarily QCed data:
for (i in 1:length(data_list)) {
	dat_name = gsub('\\..*', '', args[i])
	file_name = paste(dat_name, "qc.tsv", sep = '_')
	fwrite(data_list[[i]], file = file_name, quote = F, sep = '\t', row.names = F, col.names = T)
}

