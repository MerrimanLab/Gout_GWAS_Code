# R script to convert extremely low P values to 0 before munge_sumstats.py
#
# vroom seems to automatically convert extremely low P values into the smallest
# exponent, so this script just reads in data and writes it out

library(vroom)

args = commandArgs(trailingOnly = T)

dat = vroom(args[1])

filename = gsub('\\.txt', '.lowp.txt', args[1])

vroom_write(dat, filename)

