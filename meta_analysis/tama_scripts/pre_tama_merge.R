# Script to merge METAL output into a form ready for TAMA using MANTRA.

library(dplyr)
library(furrr)
library(purrr)
library(data.table)

# Setup multicore:
plan(strategy = 'multicore', workers = 10)
options(future.globals.maxSize= 20 * 1024 ^ 3)
options("scipen" = 999)

args = commandArgs(trailingOnly = T)

data_list = list()

dat_names = gsub('_meta.*', '', args)
dat_names = gsub('.*/', '', dat_names)

# Load data:
for (i in 1:length(args)) {
	data_list[[i]] = as.data.frame(fread(args[i]))
}

names(data_list) = dat_names

# Function to add allelic value, A1_A2, and "presence" columns to data:
add_extra = function(data) {
	minor_value = unlist(lapply(data$Allele1, function(x) switch(x, A = 1, C = 2, G = 3, T = 5, -100)))
	major_value = unlist(lapply(data$Allele2, function(x) switch(x, A = 1, C = 2, G = 3, T = 5, -100)))
	data$allelic_value = minor_value + major_value
	data$A1_A2 = paste(data$Allele1, data$Allele2, sep = '_')
	data$MarkerName = paste(data$MarkerName, data$allelic_value, sep = '_')
	rownames(data) = data$MarkerName
	return(data)
}

# Add allelic value to all data:
data_list = future_map(data_list, add_extra)

# Align all markers based on the alleles of first data set, and then merge them
# together.
#
# Since markers have alellic value added to the ID, there is no need to check
# for allelic value - just check whether the alleles match up or not. If it
# doesn't match up, flip the second data.
# Once joined via `full_join()` function, markers with differing allelic values
# will show up as independent marker (and will also prevent duplicate marker ID
# due to the allelic value).
flip_merge = function(dat1, dat2, dat2_name) {
	# First, find common markers between the two data sets:
	ind = which(dat1$MarkerName %in% dat2$MarkerName)
	snp = dat1$MarkerName[ind]
	alleles = dat1$A1_A2[ind]
	tmp = dat2[snp, ]
	if (any(snp != tmp$MarkerName)) stop('Something went wrong with pulling out common markers')
	# Find markers to flip:
	flip_ind = which(alleles != tmp$A1_A2)
	# Flip markers:
	tmp_allele = tmp$Allele1[flip_ind]
	tmp$Allele1[flip_ind] = tmp$Allele2[flip_ind]
	tmp$Allele2[flip_ind] = tmp_allele
	tmp$Freq1[flip_ind] = 1 - tmp$Freq1[flip_ind]
	tmp$Effect[flip_ind] = -tmp$Effect[flip_ind]
	tmp$A1_A2[flip_ind] = paste(tmp$Allele1, tmp$Allele2, sep = '_')
	# Add the new values:
	dat2[snp, ] = tmp
	# Merge the two data sets and return:
	merged = full_join(dat1, dat2, by = c('CHR', 'POS', 'MarkerName', 'Allele1', 'Allele2', 'allelic_value', 'A1_A2'), suffix = c('', paste0('.', dat2_name)))
	return(merged)
}

full_merge  = data_list[[1]]

for (i in 2:length(data_list)) {
	full_merge = flip_merge(full_merge, data_list[[i]], dat2_name = names(data_list)[i])
}

# Instead of allelic value, use A1_A2 values instead:
full_merge$MarkerName = paste(full_merge$CHR, full_merge$POS, full_merge$A1_A2, sep = '_')

# Remove allelic value and A1_A2 columns, and rename some columns:
full_merge = full_merge[, !((colnames(full_merge) %in% c('allelic_value', 'A1_A2')))]
colnames(full_merge)[6:10] = paste(colnames(full_merge)[6:10], names(data_list)[1], sep = '.')

# Organise output file name:
sample_group = ifelse(grepl("full", args[1]), "full/", ifelse(grepl("female", args[1]), "female/", "male/"))
out_dir = 'data/tama_dat/'
out_dir = paste(out_dir, sample_group, sep = '')
out_name = gsub('\\/', '', sample_group)
out_name = paste(out_name, 'merged_dat.tsv', sep = '_')
out_name = paste(out_dir, out_name, sep = '')

# Save the merged data:
write.table(full_merge, file = out_name, quote = F, sep = '\t', row.names = F, col.names = T)
