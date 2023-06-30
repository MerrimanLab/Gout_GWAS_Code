# Script to overlap the lead variants and Activity-by-contact (ABC) data
library(dplyr)
library(purrr)
library(vroom)

args = commandArgs(trailingOnly = T)

# Load ABC data
abc = vroom('/Volumes/archive/merrimanlab/reference_files/ABC_model/data/AllPredictions.AvgHiC.ABC0.015.minus150.ForABCPaperV3.txt')

# Load data that contains SNP, CHR, and POS
loci = vroom(args[1])
colnames(loci) = toupper(colnames(loci))
loci$CHR = as.character(loci$CHR)
loci$CHR[loci$CHR == '23'] = 'X'

# Pull out overlapping ABC sites:
res_list = list()

for (i in 1:nrow(loci)) {
	query_chr = paste('chr', loci$CHR[i], sep = '')
	query_pos = loci$POS[i]
	tmp = abc %>% filter(chr == query_chr, start <= query_pos, end >= query_pos) %>% arrange(desc(ABC.Score)) %>% slice(1:20)
	if (nrow(tmp) > 0) {
		tmp$SNP = loci$SNP[i]
		tmp$SNP.chr = query_chr
		tmp$SNP.pos = query_pos
	}
	res_list[[i]] = tmp
}

res = map_dfr(res_list, ~ .x)

res = res %>% select(contains('snp'), chr:CellType)

# Save
write.table(res, args[2], sep = '\t', col.names = T, row.names = F, quote = F)
