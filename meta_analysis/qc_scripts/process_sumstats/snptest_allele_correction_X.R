#!/bin/Rscript --vanilla

# SNPTEST gives you the MAF of allele A or B - whichever the minor one is.
# However, it will NOT tell you which allele the MAF is for, so you have to
# figure it out yourself.
#
# This is a script to solve this problem.

library(data.table)

args = commandArgs(trailingOnly = T)

data = fread(args[1])

colnames(data)[40] = 'frequentist_add_beta_1'

# For chromosome X, use the sum of A/AA and B/BB to figure out which allele is
# the minor one
#
# First identify which variants to flip:
# Note that I am going to keep the variants that have MAF = 0.5 as is (i.e. no flipping).
data$flip = (data$all_A + data$all_AA) < (data$all_BB + data$all_B)
data$flip_case = (data$cases_A + data$cases_AA) < (data$cases_BB + data$cases_B)
data$flip_control = (data$controls_A + data$controls_AA) < (data$controls_BB + data$controls_B)

# "Flip" the MAF of the wrong variants:
data$maf = ifelse(data$flip, data$all_maf, 1 - data$all_maf)
data$maf_case = ifelse(data$flip, data$cases_maf, 1 - data$cases_maf)
data$maf_control = ifelse(data$flip, data$controls_maf, 1 - data$controls_maf)

# Flip the effect size:
data$effect = ifelse(data$flip, data$frequentist_add_beta_1, -data$frequentist_add_beta_1)

# Note down how many variants had MAF = 0.5:
maf_half = which((data$all_A + data$all_AA) == (data$all_B + data$all_BB))

text = paste('There were', length(maf_half))
text = paste(text, 'with MAF = 0.5')

message(text)

# Check for NAs in the effect column for all the variants:
ind = which(is.na(data$effect))

if (length(ind) > 0) {
	data = data[-ind,]
}

# Generate a column of OR, just in case:
data$OR = exp(data$effect)

# Save file:
filename = gsub('\\.txt', '', args[1])
filename = gsub('\\.out', '', filename)
filename = paste(filename, 'corrected.txt', sep = "_")

fwrite(data, filename, col.names = T, row.names = F, quote = F, sep = ' ', na = NA, nThread = 10)

