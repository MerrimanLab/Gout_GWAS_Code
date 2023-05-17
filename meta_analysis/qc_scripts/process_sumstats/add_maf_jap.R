# Short R script to add overall MAF information and sample numbers to the
# Japanese Gout summary statistics
#
# Make sure this script is run from the root directory of the project.

# Load library:
library(data.table)
library(dplyr)

# Load data:
japonica = as.data.frame(fread('data/summary/EAS/Japanese_Japonica.tsv'))
illumina = as.data.frame(fread('data/summary/EAS/Japanese_Illumina.tsv'))

# Add sample size:
japonica$N_case = 1028
japonica$N_control = 1125
japonica$N_total = 2153

illumina$N_case = 2025
illumina$N_control = 4512
illumina$N_total = 6537

# Calculate MAC first:
japonica$MAC_case = 2 * pmin(japonica$EAF_CASE, 1 - japonica$EAF_CASE) * japonica$N_case
japonica$MAC_control = 2 * pmin(japonica$EAF_CTRL, 1 - japonica$EAF_CTRL) * japonica$N_control
japonica$MAC_total = japonica$MAC_case + japonica$MAC_control

illumina$MAC_case = 2 * pmin(illumina$EAF_CASE, 1 - illumina$EAF_CASE) * illumina$N_case
illumina$MAC_control = 2 * pmin(illumina$EAF_CTRL, 1 - illumina$EAF_CTRL) * illumina$N_control
illumina$MAC_total = illumina$MAC_case + illumina$MAC_control

# Calculate MAF:
japonica$MAF_total = (japonica$MAC_total) / (2 * japonica$N_total)
illumina$MAF_total = (illumina$MAC_total) / (2 * illumina$N_total)

# Note that the MAC is always going to be below 50%, and therefore MAF will
# also be < 0.5, which may be inconsistent with case/control MAF.
# Theoretically, overall MAF should be in between case/control MAF, so check
# for MAF that aren't in between and flip them:
tmp_maf = as.matrix(japonica[,c('EAF_CASE','EAF_CTRL', 'MAF_total')])
tmp_maf = cbind(tmp_maf, abs(tmp_maf[,1] - tmp_maf[,2]))
tmp_maf = cbind(tmp_maf, abs(tmp_maf[,1] - tmp_maf[,3]))
tmp_maf = cbind(tmp_maf, abs(tmp_maf[,2] - tmp_maf[,3]))
tmp_maf = cbind(tmp_maf, abs(tmp_maf[,4] - (tmp_maf[,5] + tmp_maf[,6])))
ind = which(!near(tmp_maf[,7], 0))
japonica$MAF_total[ind] = 1 - japonica$MAF_total[ind]

tmp_maf = as.matrix(illumina[,c('EAF_CASE','EAF_CTRL', 'MAF_total')])
tmp_maf = cbind(tmp_maf, abs(tmp_maf[,1] - tmp_maf[,2]))
tmp_maf = cbind(tmp_maf, abs(tmp_maf[,1] - tmp_maf[,3]))
tmp_maf = cbind(tmp_maf, abs(tmp_maf[,2] - tmp_maf[,3]))
tmp_maf = cbind(tmp_maf, abs(tmp_maf[,4] - (tmp_maf[,5] + tmp_maf[,6])))
ind = which(!near(tmp_maf[,7], 0))
illumina$MAF_total[ind] = 1 - illumina$MAF_total[ind]

# I used a very convoluted way to find variants to flip. Here's an alternative
# way to do it, but this method can't/doesn't handle variants when
# case/control/overall MAFs are equal to one another, so I guess the above
# method is slightly more robust.
# (Although these variants will likely be non-significant in the GWAS anyway,
# so you could probably filter these out)
#
# ind = which(!(japonica$MAF_total > pmin(japonica$EAF_CASE, japonica$EAF_CTRL) & japonica$MAF_total < pmax(japonica$EAF_CASE, japonica$EAF_CTRL)))
# japonica$MAF_total[ind] = 1 - japonica$MAF_total[ind]

# Save data:
fwrite(japonica, 'data/summary/EAS/jap_japonica_maf.tsv', quote = F, sep = '\t', row.names = F, col.names = T)
fwrite(illumina, 'data/summary/EAS/jap_illumina_maf.tsv', quote = F, sep = '\t', row.names = F, col.names = T)
