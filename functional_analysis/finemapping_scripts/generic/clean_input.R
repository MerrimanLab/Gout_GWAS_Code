# Script to pull out the correct variants for the locus and flip the direction
# of the summary stats if necessary

library(dplyr)
library(vroom)
library(purrr)

args = commandArgs(trailingOnly = T)

ld_mat = read.table(args[1], header = F, stringsAsFactors = F) %>% as.matrix

bim = read.table(args[2], header = F, stringsAsFactors = F)
colnames(bim) = c('CHR', 'SNP', 'cM', 'BP', 'A1', 'A2')

freq = read.table(args[3], header = T, stringsAsFactors = F)

z_dat = vroom(args[4])

# Remove variants with missing LD:

# First remove the variants with missing LD with every other variants:
ind = which(apply(ld_mat, 1, function(x) all(is.na(x))))

if (length(ind) > 0) {
	ld_mat = ld_mat[-ind, -ind]
	bim = bim[-ind, ]
	freq = freq[-ind, ]
}

# Remove multi-allelic variants
ind = c(which(nchar(bim$A1) > 1), which(nchar(bim$A2) > 1)) %>% unique

if (length(ind) > 0) {
	ld_mat = ld_mat[-ind, -ind]
	bim = bim[-ind, ]
	freq = freq[-ind, ]
}

# For each row, count up how many missing data there are, then remove the most
# missing variants first until there are no missing data.
# This is done to prevent unnecessary removal of variants (especially lead
# variants) that only have a single missing data point
ind = apply(ld_mat, 1, function(x) which(is.na(x)) %>% length)

while (max(ind) > 0) {
	rm = which(ind == max(ind))
	ld_mat = ld_mat[-rm, -rm]
	bim = bim[-rm, ]
	freq = bim[-rm, ]
	ind = apply(ld_mat, 1, function(x) which(is.na(x)) %>% length)
}

# Create unique ID:
bim$minor_val = unlist(lapply(bim$A1, function(x) switch(x, 'A' = 1, 'C' = 2, 'G' = 3, 'T' = 5, 100)))
bim$major_val = unlist(lapply(bim$A2, function(x) switch(x, 'A' = 1, 'C' = 2, 'G' = 3, 'T' = 5, 100)))
bim$av = bim$minor_val + bim$major_val
bim$cpid = paste(bim$CHR, bim$BP, sep = '_')
bim$avid = paste(bim$cpid, bim$av, sep = '_')

z_dat$minor_val = unlist(lapply(z_dat$minor, function(x) switch(x, 'A' = 1, 'C' = 2, 'G' = 3, 'T' = 5, 100)))
z_dat$major_val = unlist(lapply(z_dat$major, function(x) switch(x, 'A' = 1, 'C' = 2, 'G' = 3, 'T' = 5, 100)))
z_dat$av = z_dat$minor_val + z_dat$major_val
z_dat$cpid = paste(z_dat$CHR, z_dat$BP, sep = '_')
z_dat$avid = paste(z_dat$cpid, z_dat$av, sep = '_')

# Make sure both LD and bim has the same variants as the summary stats:
ind = which(!(bim$avid %in% z_dat$avid))

if (length(ind) > 0) {
	bim = bim[-ind,]
	freq = freq[-ind,]
	ld_mat = ld_mat[-ind, -ind]
}

z_dat = z_dat[z_dat$avid %in% bim$avid, ]

# Check if correct variants have been pulled out
if (length(unique(c(nrow(bim), nrow(ld_mat), nrow(z_dat)))) != 1) {
	stop('Something went wrong with getting common variants')
}

# Check order of SNPs
if (any(z_dat$avid != bim$avid)) {
	stop('Z-score data and bim file is in different order')
}

# Flip Z-scores if the alleles are different:
#
# A1 is "usually" the minor allele in PLINK, so compare this allele with the
# effect allele (minor allele) in the Z-score data
flip = which(bim$A1 != z_dat$minor)

if (length(flip) > 0) {
	tmp_allele = z_dat$minor[flip]
	z_dat$minor[flip] = z_dat$major[flip]
	z_dat$major[flip] = tmp_allele
	z_dat$effect[flip] = -z_dat$effect[flip]
	z_dat$MAF[flip] = 1 - z_dat$MAF[flip]
	z_dat$z[flip] = -z_dat$z[flip]
}

# If the summary MAF is too different to the reference MAF (greater than 5%
# discrepancy), remove from the data
maf_diff = z_dat$MAF - freq$MAF
ind = which(abs(maf_diff) > 0.05)

if (length(ind) > 0) {
	ld_mat = ld_mat[-ind, -ind]
	z_dat = z_dat[-ind,]
}

# Select only relevant columns and rename them:
program = args[5]

# For PAINTOR, generate extra table for making annotation of the locus.
# Also, make a data set with z-scores that use "normalised" beta/effect
if (toupper(program) == 'PAINTOR') {
	z_dat = z_dat %>% mutate(std_z = (z - mean(z)) / sd(z))
	z_dat = z_dat %>% select(CHR:major, effect, z, std_z)
	for_annot = z_dat %>% mutate(CHR = as.character(CHR), CHR = ifelse(CHR == '23', 'X', CHR), CHR = paste('chr', CHR, sep = '')) %>% select(CHR, BP)
	write.table(for_annot, gsub('txt', 'for_annot.txt', args[4]), sep = ' ', col.names = T, row.names = F, quote = F)
} else if (toupper(program) == 'FINEMAP') {
	# Calculate median P-value of the significant variants
	med_p = median(z_dat$P[z_dat$P <= 1e-6])
	mean_p = mean(z_dat$P[z_dat$P <= 1e-6])
	writeLines(as.character(c(med_p, mean_p)), gsub('txt', 'p_info', args[4]))
	z_dat = z_dat %>% select(SNP, CHR, BP, minor, major, MAF, effect, SE)
	colnames(z_dat) = c("rsid", "chromosome", "position", "allele1", "allele2", "maf", "beta", "se")
	# Need to correct for variants with MAF > 0.5 (i.e. those that are not
	# "minor") for FINEMAP to run properly. Do this by flipping both MAF and
	# effect/beta, and add the flip column:
	ind = which(z_dat$maf > 0.5)
	z_dat$flip = 0
	if (length(ind) > 0) {
		z_dat$beta[ind] = -z_dat$beta[ind]
		z_dat$maf[ind] = 1 - z_dat$maf[ind]
		z_dat$flip[ind] = 1
	}
}

# Save file:
write.table(ld_mat, gsub('txt', 'ld', args[4]), sep = ' ', col.names = F, row.names = F, quote = F)
write.table(z_dat, gsub('txt', 'clean.txt', args[4]), sep = ' ', col.names = T, row.names = F, quote = F)
