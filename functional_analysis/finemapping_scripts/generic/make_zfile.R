# Script to pull out region from the summary stats and generate Z scores

library(vroom)
library(purrr)
library(tidyr)
library(dplyr)

args = commandArgs(trailingOnly = T)

dat = vroom(args[1])

loci = read.table(args[2], sep ='\t', header = T, stringsAsFactors = F)
loci = loci %>% filter(EUR == 'EUR')

sex = args[3]

# For each locus, determine the locus boundary for fine-mapping. If the region
# from the boundary is >1Mb, then use 1Mb around the lead variant so that it
# doesn't extend beyond the locus boundary.

loci$locus = gsub('MB.*', '', loci$locus)
loci$locus = gsub('chr', '', loci$locus)

loci = loci %>% separate(locus, into = c('chr', 'start', 'end'), sep = '_')
loci$chr = as.numeric(loci$chr)
loci$start = as.numeric(loci$start) * 1000000
loci$end = as.numeric(loci$end) * 1000000

loci$size = loci$end - loci$start

ind = which(loci$size >= 1000000)

for (i in 1:length(ind)) {
	tmp_start = loci$BP[ind[i]] - 500000
	tmp_end = loci$BP[ind[i]] + 500000
	if (tmp_start < loci$start[ind[i]]) {
		loci$end[ind[i]] = loci$start[ind[i]] + 1000000
	} else if (tmp_end > loci$end[ind[i]]) {
		loci$start[ind[i]] = loci$end[ind[i]] - 1000000
	} else {
		loci$start[ind[i]] = tmp_start
		loci$end[ind[i]] = tmp_end
	}
}

# Save loci regions for LD
options(scipen = 100)
loci$ld_region = ifelse(loci$chr < 10, paste('0', loci$chr, sep = ''), loci$chr)
loci$ld_region = paste(loci$ld_region, as.character(loci$start), sep = ':')
loci$ld_region = paste(loci$ld_region, as.character(loci$end), sep = '-')
ld_region = loci %>% select(SNP, ld_region)

file_name = paste('data/ukbb_ld', sex, 'loci_regions.txt', sep = '/')
write.table(ld_region, file_name, sep = '\t', col.names = F, row.names = F, quote = F)

# Function to pull out given region
pull_region = function(dat, chr, start, end) {
    res = dat %>% filter(CHR == chr, between(BP, start, end))
    return(res)
}

regions = pmap(list(loci$chr, loci$start, loci$end), ~ pull_region(dat, ..1, ..2, ..3))

# Function to calculate Z-score
calc_z = function(dat) {
    dat$z = dat$effect / dat$SE
    return(dat)
}

z_dat = map(regions, ~ calc_z(.x))

z_dat = map(z_dat, ~ .x %>% select(SNP, CHR, BP, minor, major, MAF, effect, SE, z, P))

# Save Z-scores into separate files
variants = loci$SNP
out_name = paste('data/fine_mapping/finemap_dat', sex, variants, sep = '/')
out_name = paste(out_name, '.txt', sep = '')

imap(z_dat, ~ write.table(.x, out_name[.y], col.names = T, row.names = F, sep = '\t', quote = F))

out_name = gsub('finemap', 'paintor', out_name)

imap(z_dat, ~ write.table(.x, out_name[.y], col.names = T, row.names = F, sep = '\t', quote = F))

