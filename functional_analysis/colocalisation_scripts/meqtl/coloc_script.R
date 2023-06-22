# Run coloc for GWAS/meQTL
#
# Usage:
# Rscript src/coloc_script.R <path_to>/gwas_results.txt <path_to>/mqtl_results.txt

library(dplyr)
library(purrr)
library(stringr)
library(vroom)
library(coloc)

source("/Volumes/scratch/merrimanlab/riku/merrimanlab_scratch/LocusZooms/functions/locus_zoom.R")
gene_dat = read.table('/Volumes/scratch/merrimanlab/riku/merrimanlab_scratch/LocusZooms/Gencode_GRCh37_Genes_UniqueList2021.txt', sep = '\t', header = T, stringsAsFactors = F)

args = commandArgs(trailingOnly = T)

gwas = vroom(args[1], col_types = cols(P = col_character()))
mqtl = vroom(args[2], col_types = cols(pval = col_character()))

colnames(gwas) = c("CHR", "BP", "SNP", "minor", "major", "MAF", "effect", "SE", "P", "N", "cpid")
colnames(mqtl) = c("cpg", "CHR", "BP", "cpid", "effect", "SE", "P", "N", "minor", "major", "MAF", "freq_se")

gwas$varbeta = gwas$SE ^ 2
mqtl$varbeta = mqtl$SE ^ 2

# Make allelic value ID
get_av = function(a1, a2) {
	val1 = unlist(lapply(a1, function(x) switch(x, 'A' = 1, 'C' = 2, 'G' = 3, 'T' = 5, -100)))
	val2 = unlist(lapply(a2, function(x) switch(x, 'A' = 1, 'C' = 2, 'G' = 3, 'T' = 5, -100)))
	av = val1 + val2
	return(av)
}

gwas$av = get_av(gwas$minor, gwas$major)
gwas$avid = paste(gwas$cpid, gwas$av, sep = '_')

mqtl$av = get_av(mqtl$minor, mqtl$major)
mqtl$avid = paste(mqtl$cpid, mqtl$av, sep = '_')

# Pull out common variants
common = gwas$avid[which(gwas$avid %in% mqtl$avid)]

gwas = gwas[which(gwas$avid %in% common),] %>% arrange(CHR, BP)
mqtl = mqtl[which(mqtl$avid %in% common),] %>% arrange(CHR, BP)

# Make sure all the variants are +/-500kb from the GWAS lead SNP, but first
# need to figure out what the lead SNP is.
# The lead SNP is in the filename, so find the cpid in the filename
ind = map_lgl(gwas$cpid, ~ grepl(.x, args[1]))
ind = which(ind)

snp = gwas$SNP[ind]
cpg = mqtl$cpg[ind]
chr = gwas$CHR[ind]
pos = gwas$BP[ind]
start = pos - 500000
end = pos + 500000

# Now filter out the variants outside the +/-500kb window
gwas = gwas %>% filter(between(BP, start, end))
mqtl = mqtl %>% filter(between(BP, start, end))

if (nrow(gwas) != nrow(mqtl)) {
	stop("Something went wrong with merging data sets - ", gsub('.gwas.txt', '', args[1]))
}

# Align alleles
ind = which(gwas$minor != mqtl$minor)

if (length(ind) > 0) {
	tmp = mqtl$minor[ind]
	mqtl$minor[ind] = mqtl$major[ind]
	mqtl$major[ind] = tmp
	mqtl$MAF[ind] = 1 - mqtl$MAF[ind]
	mqtl$effect[ind] = -mqtl$effect[ind]
}

# (Approximate) Proportion of cases in the EUR meta-analysis data
prop = 100661 / 2106003

# Make coloc input lists
gwas_list = list(snp = gwas$SNP, beta = gwas$effect, MAF = gwas$MAF, N = gwas$N, varbeta = gwas$varbeta, type = "cc", s = prop)
mqtl_list = list(snp = gwas$SNP, beta = mqtl$effect, MAF = mqtl$MAF, N = mqtl$N, varbeta = mqtl$varbeta, type = "quant")

res = coloc.abf(gwas_list, mqtl_list)

# Make an informative rowname with rsID and cpgID

# Load in 450K CHR/POS information
cpg_annot = vroom('data/meqtl_data/cpg450k_chrpos.txt')

cpg.chr = cpg_annot$CHR[which(cpg_annot$cpg == cpg)]
cpg.pos = cpg_annot$BP[which(cpg_annot$cpg == cpg)]

sum_res = t(as.data.frame(c(SNP = snp, SNP.CHR = chr, SNP.BP = pos, mqtl_cpg = cpg, cpg.CHR = cpg.chr, cpg.BP = cpg.pos, res$summary)))

# Save results and combined summary stats
res_name = gsub('gwas.txt', 'coloc_summary.txt', args[1])
write.table(sum_res, res_name, sep = '\t', col.names = T, row.names = F, quote = F)

sumstat = left_join(gwas, mqtl, by = c('CHR', 'BP', 'cpid', 'minor', 'major', 'av', 'avid'), suffix = c('.gwas', '.mqtl'))

summary_name = gsub('gwas.txt', 'merged_summary.txt', args[1])
write.table(sumstat, summary_name, sep = '\t', col.names = T, row.names = F, quote = F)
