################################################################################
# Rscript for reporting the number of variants violating defined thresholds
# and outputs a file containing the violating variants
################################################################################

library(vroom)
library(dplyr)

source('src/across_study_qc.func.R')

args = commandArgs(trailingOnly = T)

data_list = list()

dat_names = gsub('_clean.*', '', args)
dat_names = gsub('.*/', '', dat_names)

# Load data:
for (i in 1:length(args)) {
	data_list[[i]] = vroom(args[i])
}

names(data_list) = dat_names

reorder = order(unlist(lapply(data_list, nrow)), decreasing = T)

data_list = data_list[reorder]
dat_names = dat_names[reorder]

# Remove and report SNPs that are bad (in terms of allele coding) between data
# sets:
data_list = remove_bad_alleles(data_list)

################################################################################
# Code to fix the "cross" in MAF vs. MAF plot
#
# Since different set of variants may overlap between different studies, it is
# important to compare the frequencies for all the data sets, not just with the
# reference/first data set.
#
# If these are not checked over all the data sets, there is no guarantee that
# the meta-analysis result is reliable in the case of the variants that are
# not present in the reference data set, but present in other data sets.
# E.g. study1 (ref) might not have the variant, but study2 and study3 do, and
# they have conflicting MAF/beta -> result may be misleading as MAF and beta is
# for the wrong allele in one of the study.
#
# This will be fixed by iteratively flipping the alleles based on the reference
# data set, then the next data set in the list, and so on. Since only the
# comparison data set is flipped, there will be no "re-flipping" occuring, as
# long as those variants that have already been flipped are removed
# subsequently.

# Generate all the combinations to check the alleles against:
comb = combn(length(data_list), 2)

clean_dat = data_list

ign = NULL

# Correct the alleles:
for (i in 1:ncol(comb)) {
	dat1 = clean_dat[[comb[1,i]]]
	dat2 = clean_dat[[comb[2,i]]]
	data_list[[comb[2,i]]] = standardise_allele(dat1, dat2, tol_avg = 0.01, tol_diff = 0.01, ignore = ign)
	if (comb[2,i] == length(clean_dat)) {
		ign = unique(c(ign, clean_dat[[comb[1,i]]]$cpid))
	}
}

################################################################################
# Remove the extreme variants from data, using the list generated earlier
################################################################################

# Generate vector of output name:
ancestry = gsub('\\/[^\\/]*$', '', args)
ancestry = unique(gsub('.*/', '', ancestry))

sample_group = ifelse(grepl("full", args[1]), "full/", ifelse(grepl("female", args[1]), "female/", "male/"))

rds_dir = paste('results/pre_meta_qc', ancestry, sep = '/')
rds_dir = paste(rds_dir, sample_group, sep = '/')

rds_file = paste(rds_dir, dat_names, sep = '')
rds_file = paste(rds_file, '_extreme_snp_list.rds', sep = '')

# Remove MAF < 0.01% and absolute effect size > 10:
for (i in 1:length(clean_dat)) {
	tmp_dat = readRDS(rds_file[i])
	if (is.vector(tmp_dat)) {
		next()
	}
	ex_snps = tmp_dat[which(abs(tmp_dat$effect) > 10 | tmp_dat$MAF < 0.001 | tmp_dat$MAF > 0.999), ]
	ex_snps = unique(paste(ex_snps$cpid, ex_snps$allelic_value, sep = '_'))
	tmp_res = clean_dat[[i]]
	tmp_res = tmp_res[which(!(tmp_res$cpavid %in% ex_snps)),]
	clean_dat[[i]] = tmp_res[which(!(tmp_res$cpavid %in% ex_snps)),]
}

################################################################################
# Save data
################################################################################

tmp = gsub('\\..*', '', args[reorder])
file_name = paste(tmp, "final.txt", sep = '_')

for (i in 1:length(clean_dat)) {
	message('Writing ', names(clean_dat[i]), ' as ', file_name[i])
	vroom_write(clean_dat[[i]], file_name[i])
}

