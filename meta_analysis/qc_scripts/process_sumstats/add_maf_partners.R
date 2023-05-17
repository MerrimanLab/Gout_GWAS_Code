# R script to clean and calculate case/control/overall MAF of Partners data

library(dplyr)
library(readr)
library(data.table)

args = commandArgs(trailingOnly = T)

dat = read_delim(args[1], delim = '\t')

colnames(dat) = c("CHR", "POS", "SNP", "minor", "major", "MAF", "MAF_SE", "MinFreq", "MaxFreq", "effect", "SE", "P", "Direction", "HetISq", "HetChiSq", "HetDf", "HetPVal", "Total_N_REF_Case", "Total_N_REF_Control", "Total_N_HET_Case", "Total_N_HET_Control", "Total_N_ALT_Case", "Total_N_ALT_Control")

# Uppercase the alleles, calculate MACs and MAFs for case/control/total:
dat = dat %>% mutate(
			   minor = toupper(minor),
			   major = toupper(major),
			   MAC_case = (2 * Total_N_REF_Case) + Total_N_HET_Case,
			   MAC_control = (2 * Total_N_REF_Control) + Total_N_HET_Control,
			   MAC = MAC_case + MAC_control,
			   N_case = Total_N_REF_Case + Total_N_HET_Case + Total_N_ALT_Case,
			   N_control = Total_N_REF_Control + Total_N_HET_Control + Total_N_ALT_Control,
			   N = N_case + N_control,
			   MAF_case = MAC_case / (2 * N_case),
			   MAF_control = MAC_control / (2 * N_control),
			   MAF_calc = MAC / (2 * N)
			   ) %>% select(CHR:MAF, effect:P, MAC_case:MAF_calc)

# Flip case/control MAF if calculated MAF < 0.5 and true MAF > 0.5 (and vice
# versa)
dat = dat %>% mutate(
			   avg_MAF = (MAF + MAF_calc) / 2,
			   flip = near(0.5, avg_MAF, 0.05) | (MAF <= 0.475 & MAF_calc > 0.525) | (MAF >= 0.525 & MAF_calc < 0.475),
			   MAC_case = if_else(flip, (2 * N_case) - MAC_case, MAC_case),
			   MAC_control = if_else(flip, (2 * N_control) - MAC_control, MAC_control),
			   MAC = if_else(flip, (2 * N) - MAC, MAC),
			   MAF_calc = if_else(flip, 1 - MAF_calc, MAF_calc),
			   MAF_case = if_else(flip, 1 - MAF_case, MAF_case),
			   MAF_control = if_else(flip, 1 - MAF_control, MAF_control)
			   ) %>% select(CHR:MAF_control)

# Save data
filename = gsub('\\..*', '_maf.tsv', args[1])
write_tsv(dat, filename)
