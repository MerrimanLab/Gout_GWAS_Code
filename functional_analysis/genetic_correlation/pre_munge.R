# Add rsID to the summary stats and make sure the effect allele is correct
library(dplyr)
library(vroom)

args = commandArgs(trailingOnly = T)

dat = vroom(args[1], na = c('', 'NA', 'NaN'))
rsid = vroom('data/neale_ukbb/sumstats_variants.clean.txt')

# Add rsIDs and filter out any SNPs with no rsID, low confident SNPs, and no
# beta value
res = left_join(dat, rsid, by = c('variant' = 'ID')) %>% filter(!is.na(SNP)) %>% filter(!low_confidence_variant) %>% filter(!is.na(beta))

# Make sure the A1 allele is the effect allele
ind = which(res$minor_allele != res$A1)

tmp = res$A1[ind]
res$A1[ind] = res$A2[ind]
res$A2[ind] = tmp

# Double-check that the minor/effect allele is the A1 allele
if(any(res$minor_allele != res$A1)) {
	stop('At least one A1 was not a minor/effect allele')
}

res = res %>% select(SNP, CHR, POS, A1, A2, minor_AF, low_confidence_variant:pval)
colnames(res) = c("SNP", "CHR", "POS", "A1", "A2", "MAF", "low_confidence_variant", "n_samples", "AC", "ytx", "beta", "SE", "tstat", "P")

# Save:
file_name = gsub('tsv.gz', 'clean.tsv.gz', args[1])
vroom_write(res, file_name)
