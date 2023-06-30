# Check LD of clump lead variant across loci
library(dplyr)
library(vroom)
library(tidyr)
library(GenomicRanges)

ld_tab = read.delim('res/urate_clumping/clump_lead_ld.ld', sep = '', header = T, stringsAsFactors = F)
ld_tab = ld_tab %>% filter(SNP_A != SNP_B, CHR_A == CHR_B)
ld_tab = ld_tab %>% filter(R2 > 0.2)

write.table(ld_tab, 'res/urate_clumping/clump_lead_ld.r2_0.2.txt', sep = '\t', col.names = T, row.names = F, quote = F)

# Load in the preliminary clump result and merge the loci that have r2 > 0.2
crude = read.table('dat/urate_clumping/crude_clump.txt', sep = '\t', header = T, stringsAsFactors = F)
clump = read.table('res/urate_clumping/all_clump_lead.txt', sep = '\t', header = T, stringsAsFactors = F)
clump = left_join(crude, clump, by = c('simplified_name' = 'loci'))

# Add P-value to the clump lead to make it easier to track what the lead
# variant would be in the merged locus
dat = vroom('res/urate_meta/tin_ukbb.urate_meta1.for_clumping.txt')
dat = dat %>% filter(P <= 5e-8)

clump = dat %>% select(SNP, P) %>% left_join(clump, ., by = c('clump_lead' = 'SNP'))
colnames(clump)[ncol(clump)] = 'clump_lead_p'

# Merge loci that have r2 > 0.2 with a lead variant in neighbouring loci
merge_snp = unique(c(ld_tab$SNP_A, ld_tab$SNP_B))

for (i in 2:nrow(clump)){
	snp1 = clump$clump_lead[i - 1]
	snp2 = clump$clump_lead[i]
	if (snp1 %in% merge_snp && snp2 %in% merge_snp) {
		ind = c(i - 1, i)
		clump$start[ind] = min(clump$start[ind])
		clump$end[ind] = max(clump$end[ind])
		# There must be a better way, but I'm indexing to narrow down the lead
		# variants first so the which() indexes the correct variant
		clump$clump_lead[ind] = clump$clump_lead[ind][which(clump$clump_lead_p[ind] == min(clump$clump_lead_p[ind]))]
		clump$clump_lead_p[ind] = min(clump$clump_lead_p[ind])
	}
	# snps = ld_tab %>% filter(SNP_A == merge_snp[i] | SNP_B == merge_snp[i])
	 # snps = unique(c(snps$SNP_A, snps$SNP_B))
	 # ind = which(clump$clump_lead %in% snps)
	 # clump$start[ind] = min(clump$start[ind])
	 # clump$end[ind] = max(clump$end[ind])
}

# Merge loci together:
to_merge = clump %>% select(seqnames:strand) %>% GRanges %>% reduce %>% as.data.frame
to_merge = to_merge %>% mutate(loci = paste('chr', seqnames, ':', start, '-', end, sep = ''), simplified_name = gsub('[:-]', '_', loci), width_kb = ceiling(width / 1000), seqnames = as.numeric(seqnames))
write.table(to_merge, 'dat/urate_clumping/all_loci.merged_r2_0.2.txt', sep = '\t', col.names = T, row.names = F, quote = F)

# Pull out the loci that got merged:
merged = clump %>% distinct(simplified_name, clump_lead) %>% left_join(to_merge, .) %>% filter(is.na(clump_lead))
write.table(merged, 'dat/urate_clumping/reclump.txt', sep = '\t', col.names = T, row.names = F, quote = F)

# Save the loci that needs to be re-clumped
for (i in 1:nrow(merged)) {
	tmp_dat = dat %>% filter(CHR == merged$seqnames[i], between(POS, merged$start[i], merged$end[i]))
	out_file = paste('dat/urate_clumping/', merged$simplified_name[i], '.txt', sep = '')
	write.table(tmp_dat, out_file, sep = '\t', col.names = T, row.names = F, quote = F)
}

