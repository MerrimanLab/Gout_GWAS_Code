library(dplyr)
library(vroom)
library(GenomicRanges)

# Make loci based on Ruth/Tanya's method

dat = vroom('res/urate_meta/tin_ukbb.urate_meta1.for_clumping.txt')

# Filter for P <= 1e-7 and pad +/-50kb either side
dat_nom = dat %>% filter(P <= 1e-7) %>% select(SNP, CHR, POS, P, P_orig, logP) %>% mutate(start = POS - 50000, end = POS + 50000)

loci_ranges = GRanges(dat_nom) %>% reduce(.) %>% as.data.frame %>% mutate(loci = paste('chr', seqnames, ':', start, '-', end, sep = ''))

# If loci has <= 1 SNP with P <= 5e-8, discard
rm_loci = c()
for (i in 1:nrow(loci_ranges)) {
	tmp_dat = dat_nom %>% filter(CHR == loci_ranges$seqnames[i], between(POS, loci_ranges$start[i], loci_ranges$end[i]))
	count = length(which(tmp_dat$P <= 5e-8))
	if (count <= 1) {
		rm_loci = c(rm_loci, i)
	}
}

loci_ranges = loci_ranges[-(rm_loci), ]

# Save each locus into separate files
dat_clump = dat %>% filter(P <= 5e-8)

for (i in 1:nrow(loci_ranges)) {
	tmp_dat = dat_clump %>% filter(CHR == loci_ranges$seqnames[i], between(POS, loci_ranges$start[i], loci_ranges$end[i]))
	out_file = paste('dat/urate_clumping/', gsub('[:-]', '_', loci_ranges$loci[i]), '.txt', sep = '')
	write.table(tmp_dat, out_file, sep = '\t', col.names = T, row.names = F, quote = F)
}

# Add columns to simplify subsequent commands
loci_ranges$simplified_name = gsub('[:-]', '_', loci_ranges$loci)
loci_ranges$width_kb = ceiling(loci_ranges$width / 1000)

write.table(loci_ranges, 'dat/urate_clumping/crude_clump.txt', sep = '\t', col.names = T, row.names = F, quote = F)
