# Rscript that has all the functions for meta-analysis QC

# Function to fix alleles that have the wrong reference allele based on
# a reference data set, where dat1 = reference and dat2 = comparison.
#
# Only the alleles in the comparison data set, not the reference data set, will
# be altered and returned.
#
# This fixes the variants that causes the characteristic "cross" in the MAF vs
# MAF plot by flipping the alleles for the variants that have the sum of MAF
# (between the two data sets being compared) close to 1 (or the average is
# close to 0.5), since those variants in the "cross" will have inverted MAF
# (e.g. 0.3 in study1 vs 0.7 in study2).
#
# Also, for those variants that are already close to 0.5 (i.e. in the middle of
# the cross) will be excluded for the flipping, since they are already in the
# grey-zone in terms of which data set has the correct reference allele.
#
# The tol_avg adjusts the tolerance level of how close the average MAF should
# be to 0.5, and tol_diff adjusts the tolerance of how close the difference in MAF
# between the reference and comparison data sets should be to 0.
#
# The ignore option lets you to ignore certain variants, so you don't have to
# compare variants that have already been compared.
standardise_allele = function(dat1, dat2, tol_avg = 0.01, tol_diff = 0.01, ignore = NULL) {
	# Get common variants between the two studies:
	id1 = dat1$cpid
	id2 = dat2$cpid
	common_snps = unique(id1[which(id1 %in% id2)])

	# Ignore variants in ignore:
	if (!is.null(ignore)) {
		common_snps = common_snps[which(!(common_snps %in% ignore))]
	}

	# Subset the data:
	tmp1 = dat1[which(dat1$cpid %in% common_snps),]
	tmp2 = dat2[which(dat2$cpid %in% common_snps),]

	# Order the data:
	tmp1 = tmp1[order(tmp1$cpid),]
	tmp2 = tmp2[order(tmp2$cpid),]

	# Generate A1_A2 column:
	tmp1$A1_A2 = paste(tmp1$minor, tmp1$major, sep = '_')
	tmp2$A1_A2 = paste(tmp2$minor, tmp2$major, sep = '_')

	# Switch MAF, effects, and alleles based on tmp1, so they are consistent
	# between the studies:
	ind = which(tmp1$A1_A2 != tmp2$A1_A2)
	tmp2[ind,] = flip_allele(tmp2[ind,], flip_allele = T)

	# Find variants that have average MAF of 0.5 that is not in the center of
	# the cross:
	avg_maf = (tmp1$MAF + tmp2$MAF) / 2
	maf_diff = abs(tmp1$MAF - tmp2$MAF)
	ind = which(near(0.5, avg_maf, tol = tol_avg) & !near(0, maf_diff, tol = tol_diff))

	# Change the MAF/beta so it is referring to the correct allele.
	# Note that the change is made in dat2, not tmp2, since tmp2 is only
	# a subset that is common with the reference data.
	if (length(ind) > 0) {
		flip = which(dat2$cpid %in% tmp2$cpid[ind])
		dat2[flip,] = flip_allele(dat2[flip,], flip_allele = F)
	}

	# Report:
	message(length(ind), " variants were flipped due to MAF difference")

	# Return the comparison tmpa:
	return(dat2)
}

# Quick function to flip the MAF/beta for a given data set.
#
# flip_allele option allows you to decide whether to flip the minor/major
# alleles as well as MAF/beta.
#
# flip_allele should only be False if you want to change the allele the
# MAF/beta is referring to (e.g. in the case where the MAF/beta is referring to
# the wrong allele).
flip_allele = function(data, flip_allele = T) {
	if (flip_allele) {
		tmp = data$minor
		data$minor = data$major
		data$major = tmp
		data$A1_A2 = paste(data$minor, data$major, sep = '_')
	}

	data$effect = -data$effect
	data$MAF = 1 - data$MAF

	# Deal with case/control MAF, if present:
	if (any(grepl('MAF_case|MAF_control', colnames(data)))) {
		data$MAF_case = 1 - data$MAF_case
		data$MAF_control = 1 - data$MAF_control
	}

	return(data)
}

# Function to remove variants that have wrong minor/major allele pairs, based
# on the reference data set (the first data set in the input list).
#
# The variants that have the wrong alleles will be removed from the
# non-reference data set, and the variants will NOT be removed from the
# reference data set. Therefore, the reference data set should be the data set
# that you want to keep as much information about (i.e. with the largest number
# of variants and samples):
remove_bad_alleles = function(data_list) {
	for (i in 2:length(data_list)) {
		dat1 = data_list[[1]]
		dat2 = data_list[[i]]

		# Get common variants between the two studies:
		id1 = dat1$cpid
		id2 = dat2$cpid
		common_snps = unique(id1[which(id1 %in% id2)])

		# Subset the data:
		dat1 = dat1[which(dat1$cpid %in% common_snps),]
		dat2 = dat2[which(dat2$cpid %in% common_snps),]

		# Order the data:
		dat1 = dat1[order(dat1$cpid),]
		dat2 = dat2[order(dat2$cpid),]

		# Create a new ID with allelic value info:
		dat1$cpavid = paste(dat1$cpid, dat1$allelic_value, sep = '_')
		dat2$cpavid = paste(dat2$cpid, dat2$allelic_value, sep = '_')

		# Identify variants that have different allelic values (i.e. variants
		# that violate alleles):
		ind = which(!(dat1$cpavid %in% dat2$cpavid))
		bad_snps = dat1$cpavid[ind]

		# Remove bad alleles from the non-reference data:
		clean_dat = data_list[[i]]
		clean_dat$cpavid = paste(clean_dat$cpid, clean_dat$allelic_value, sep = '_')
		ind = which(!(clean_dat$cpavid %in% bad_snps))
		clean_dat = clean_dat[ind,]
		data_list[[i]] = clean_dat

		# Report:
		txt = paste("Between", names(data_list)[1], sep = ' ')
		txt = paste(txt, "and", sep = ' ')
		txt = paste(txt, names(data_list)[i], sep = ' ')
		txt = paste(txt, ", there were", sep = '')
		txt = paste(txt, length(bad_snps), sep = ' ')
		txt = paste(txt, "bad SNPs", sep = ' ')
		message(txt)
	}
	# Return final data list:
	return(data_list)
}

