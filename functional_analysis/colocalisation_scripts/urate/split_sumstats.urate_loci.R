# Extract summary stats for the urate loci from gout and urate meta-analysis
# results
library(dplyr)
library(vroom)

# Load loci info
loci = read.table('res/urate_clumping/tin_ukbb.loci_list.txt', sep = '\t', header = T, stringsAsFactors = F)
loci = loci %>% select(lead)

# Load summary stats:
gout = vroom('/Volumes/archive/merrimanlab/major_gwas_paper_archive/results_meta_analysis/EUR/full/EUR_meta_full1_clean_rsid.nfiltered.biallelic.txt')

urate = vroom('res/urate_meta/tin_ukbb.urate_meta1.clean.txt')

# Add CHR/POS info to loci
loci = urate %>% select(SNP, CHR, POS) %>% left_join(loci, ., by = c('lead' = 'SNP'))
colnames(loci)[1] = 'SNP'

# Pull out relevant regions for coloc
for (i in 1:nrow(loci)) {
	chr = loci$CHR[i]
	start = loci$POS[i] - 500000
	end = loci$POS[i] + 500000
	gout_subset = gout %>% filter(CHR == chr, between(BP, start, end))
	urate_subset = urate %>% filter(CHR == chr, between(POS, start, end))
	out_file = paste('dat/coloc_dat/TRAIT/', loci$SNP[i], '.TRAIT.txt', sep = '')
	write.table(gout_subset, gsub('TRAIT', 'gout', out_file), sep = '\t', col.names = T, row.names = F, quote = F)
	write.table(urate_subset, gsub('TRAIT', 'urate', out_file), sep = '\t', col.names = T, row.names = F, quote = F)
}
