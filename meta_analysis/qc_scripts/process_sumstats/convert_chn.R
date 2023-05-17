# Small R script to add in MAF information to the summary stats.

library(dplyr)
library(vroom)

args = commandArgs(trailingOnly = T)

dat = vroom(args[1])

# Change column names, so there's less risk of puzzling major/minor:
colnames(dat) = c("CHR", "SNP", "POS", "minor", "major", "MAF_case", "MAF_control", "OR", "SE", "P")

# Replace "<" characters (for deletion variants) in dat to "<DEL>":
ind_minor = which(dat$minor == '<')
ind_major = which(dat$major == '<')
dat[ind_minor, 'minor'] = "<DEL>"
dat[ind_major, 'major'] = "<DEL>"

# Join A1 and A2:
dat$A1_A2 = paste(dat$minor, dat$major, sep = '_')

# Add MAF info and N to summary stats (piece of code taken from
# src/add_maf_jap.R):
case = ifelse(grepl('4653', args), 4653, 4210)
control = ifelse(grepl('4653', args), 4599, 4566)
N = case + control

dat$MAC_case = 2 * pmin(dat$MAF_case, 1 - dat$MAF_case) * case
dat$MAC_control = 2 * pmin(dat$MAF_control, 1 - dat$MAF_control) * control
dat$MAC_total = dat$MAC_case + dat$MAC_control
dat$N = N
dat$N_case = case
dat$N_control = control

# Calculate MAF:
dat$MAF_total = (dat$MAC_total) / (2 * N)

# Note that the MAC is always going to be below 50%, and therefore MAF will
# also be < 0.5, which may be inconsistent with case/control MAF.
# Theoretically, overall MAF should be in between case/control MAF, so check
# for MAF that aren't in between and flip them:
tmp_maf = as.matrix(dat[,c('MAF_case','MAF_control', 'MAF_total')])
tmp_maf = cbind(tmp_maf, abs(tmp_maf[,1] - tmp_maf[,2]))
tmp_maf = cbind(tmp_maf, abs(tmp_maf[,1] - tmp_maf[,3]))
tmp_maf = cbind(tmp_maf, abs(tmp_maf[,2] - tmp_maf[,3]))
tmp_maf = cbind(tmp_maf, abs(tmp_maf[,4] - (tmp_maf[,5] + tmp_maf[,6])))
ind = which(!near(tmp_maf[,7], 0))
dat$MAF_total[ind] = 1 - dat$MAF_total[ind]

# Convert OR to beta:
dat$beta = log(dat$OR)

# Save data:
out_name = paste(args, '.converted', sep = '')

vroom_write(dat, out_name)
