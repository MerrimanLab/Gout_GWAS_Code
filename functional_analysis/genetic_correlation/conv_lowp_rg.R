# R script to convert extremely low P values to 0 before munge_sumstats.py
#
# vroom seems to automatically convert extremely low P values into the smallest
# exponent, so this script just reads in data and writes it out
#
# However, it is best to remove variants that are extremely significant for
# genetic correlation analysis, so convert any p <= 1e-300 into 0, and
# munge_sumstats.py should handle the rest.

library(vroom)

args = commandArgs(trailingOnly = T)

dat = vroom(args[1])

# Figure out which column is the P column
col_ind = which(colnames(dat) %in% c('P', 'P-value'))

# Convert extreme P to 0
ind = which(dat[, col_ind] <= 1e-300)
dat[ind, col_ind] = 0

filename = gsub('\\.txt', '.lowp.txt', args[1])

vroom_write(dat, filename)

