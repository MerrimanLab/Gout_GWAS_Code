# Filter out results based on the variants Tony is interested. These will be
# lead variants that were significant in urate, but not in gout (P > 0.01)

library(dplyr)
library(tidyr)
library(vroom)

urate = vroom('res/urate_meta/tin_ukbb.urate_meta1.clean.txt')
gout = vroom('/Volumes/archive/merrimanlab/major_gwas_paper_archive/results_meta_analysis/EUR/full/EUR_meta_full1_clean_rsid.nfiltered.biallelic.txt')

# Load coloc result for urate loci:
coloc = read.table('res/coloc/gout_urate/coloc_res.urate_loci.txt', sep = '\t', header = T, stringsAsFactors = F)

tmp_res = coloc %>% filter(PP.H4.abf < 0.8) %>% mutate(h2_plus_h4 = PP.H2.abf + PP.H4.abf) %>% filter(PP.H2.abf >= 0.8 | PP.H3.abf >= 0.8 | h2_plus_h4 >= 0.8) %>% mutate(cpid = paste(CHR, POS, sep = '_'))

# Add urate and gout P values
tmp_res = urate %>% filter(cpid %in% tmp_res$cpid) %>% select(CHR, POS, P) %>% left_join(tmp_res, ., by = c('CHR', 'POS'))
tmp_res = gout %>% filter(cpid %in% tmp_res$cpid) %>% select(CHR, BP, P) %>% left_join(tmp_res, ., by = c('CHR', 'POS' = 'BP'), suffix = c('.urate', '.gout'))

# For each locus that didn't have urate lead variant in the gout GWAS, I need
# to figure out what the lead variant was in the region after merging urate and
# gout summary stats (since only the common variants between the two were used
# for coloc)
no_p = tmp_res %>% filter(is.na(P.gout)) %>% pull(SNP)
top_snp = c()

for (i in 1:length(no_p)) {
	file_pat = 'dat/coloc_dat/TRAIT/SNP.TRAIT.txt'
	file_pat = gsub('SNP', no_p[i], file_pat)

	urate_file = gsub('TRAIT', 'urate', file_pat)
	tmp_urate = read.table(urate_file, sep = '\t', header = T, stringsAsFactors = F)
	tmp_urate = tmp_urate %>% select(CHR:freq1, effect:P, N)

	gout_file = gsub('TRAIT', 'gout', file_pat)
	tmp_gout = read.table(gout_file, sep = '\t', header = T, stringsAsFactors = F)
	tmp_gout = tmp_gout %>% select(CHR:BP, cpid, minor:MAF, effect:P, N)
	colnames(tmp_gout) = colnames(tmp_urate)

	# Get common variants:
	common = c(tmp_urate$cpid, tmp_gout$cpid)
	common = common[which(duplicated(common))]

	tmp_urate = tmp_urate %>% filter(cpid %in% common)
	tmp_gout = tmp_gout %>% filter(cpid %in% common)

	# Check the alleles are the same:
	tmp_urate$av1 = unlist(lapply(tmp_urate$allele1, function(x) switch(x, A = 1, C = 2, G = 3, T = 5)))
	tmp_urate$av2 = unlist(lapply(tmp_urate$allele2, function(x) switch(x, A = 1, C = 2, G = 3, T = 5)))
	tmp_urate$av = tmp_urate$av1 + tmp_urate$av2

	tmp_gout$av1 = unlist(lapply(tmp_gout$allele1, function(x) switch(x, A = 1, C = 2, G = 3, T = 5)))
	tmp_gout$av2 = unlist(lapply(tmp_gout$allele2, function(x) switch(x, A = 1, C = 2, G = 3, T = 5)))
	tmp_gout$av = tmp_gout$av1 + tmp_gout$av2

	ind = which(tmp_urate$av != tmp_gout$av)

	if(length(ind) > 0) {
		tmp_urate = tmp_urate[-(ind),]
		tmp_gout = tmp_gout[-(ind),]
	}

	tmp_snp = tmp_urate %>% arrange(P) %>% pull(cpid)
	top_snp = c(top_snp, tmp_snp[1])
}

# Make it in to a data frame so it's easier to merge
diff_lead = data.frame(lead_urate = no_p, lead_coloc = top_snp)

diff_lead = left_join(diff_lead, urate, by = c('lead_coloc' = 'cpid')) %>% select(lead_urate:lead_coloc, SNP, CHR:POS, P)
diff_lead = gout %>% select(SNP, P) %>% left_join(diff_lead, ., by = 'SNP', suffix = c('.urate', '.gout'))

# Merge this info to the other table:
tmp_res = diff_lead %>% select(lead_urate, SNP:P.gout) %>% left_join(tmp_res, ., by = c('SNP' = 'lead_urate'), suffix = c('', '.coloc_lead'))

# Check for the loci that passes Tony's definition
# These variants will be used to filter out the relevant loci/variants from the
# urate/GTEx cis-eQTL coloc results
res = tmp_res %>% filter(P.gout > 0.01 | P.gout.coloc_lead > 0.01)

write.table(res, 'res/coloc/gout_urate/coloc_res.h2_h3_filtered.txt', sep = '\t', col.names = T, row.names = F, quote = F)

