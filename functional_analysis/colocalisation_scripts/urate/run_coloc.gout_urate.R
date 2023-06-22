# Run coloc of gout and urate, given the SNP of interest
library(dplyr)
library(tidyr)
library(coloc)

args = commandArgs(trailingOnly = TRUE)

snp = args[1]

file_pat = 'dat/coloc_dat/TRAIT/SNP.TRAIT.txt'
file_pat = gsub('SNP', snp, file_pat)

urate_file = gsub('TRAIT', 'urate', file_pat)
urate = read.table(urate_file, sep = '\t', header = T, stringsAsFactors = F)
urate = urate %>% select(CHR:freq1, effect:P, N)

gout_file = gsub('TRAIT', 'gout', file_pat)
gout = read.table(gout_file, sep = '\t', header = T, stringsAsFactors = F)
gout = gout %>% select(CHR:BP, cpid, minor:MAF, effect:P, N)
colnames(gout) = colnames(urate)

# Get common variants:
common = c(urate$cpid, gout$cpid)
common = common[which(duplicated(common))]

urate = urate %>% filter(cpid %in% common)
gout = gout %>% filter(cpid %in% common)

# Check the alleles are the same:
urate$av1 = unlist(lapply(urate$allele1, function(x) switch(x, A = 1, C = 2, G = 3, T = 5)))
urate$av2 = unlist(lapply(urate$allele2, function(x) switch(x, A = 1, C = 2, G = 3, T = 5)))
urate$av = urate$av1 + urate$av2

gout$av1 = unlist(lapply(gout$allele1, function(x) switch(x, A = 1, C = 2, G = 3, T = 5)))
gout$av2 = unlist(lapply(gout$allele2, function(x) switch(x, A = 1, C = 2, G = 3, T = 5)))
gout$av = gout$av1 + gout$av2

ind = which(urate$av != gout$av)

if(length(ind) > 0) {
	urate = urate[-(ind),]
	gout = gout[-(ind),]
}

# Align the alleles
ind = which(urate$allele1 != gout$allele1)

if(length(ind) > 0) {
	tmp = gout$allele1[ind]
	gout$allele1[ind] = gout$allele2[ind]
	gout$allele2[ind] = tmp
	gout$effect[ind] = -(gout$effect[ind])
	gout$freq1[ind] = 1 - gout$freq1[ind]
}

# Function to pull out relevant coloc info from summary stats and run coloc
# Case proportion for EUR full GWAS = 100661 / 2106003 = 0.0478
run_coloc = function(gout, urate, prop_case = 0.0478) {
	gout = list(snp = gout$cpid,
				beta = gout$effect,
				MAF = gout$freq1,
				N = gout$N,
				varbeta = (gout$SE ^ 2),
				pvalues = gout$P,
				type = 'cc',
				s = prop_case)

	urate = list(snp = urate$cpid,
			   beta = urate$effect,
			   MAF = urate$freq1,
			   N = urate$N,
			   varbeta = (urate$SE ^ 2),
			   pvalues = urate$P,
			   type = 'quant')

	results = coloc.abf(gout, urate)
	return(results$summary)
}

coloc_res = run_coloc(gout, urate)
coloc_res = as.data.frame(t(coloc_res))
coloc_res$SNP = snp
coloc_res = coloc_res %>% select(SNP, nsnps:PP.H4.abf)

# Save result
out_name = paste('res/coloc/gout_urate/', snp, '.coloc_res.txt', sep = '')
write.table(coloc_res, out_name, sep = '\t', col.names = T, row.names = F, quote = F)
