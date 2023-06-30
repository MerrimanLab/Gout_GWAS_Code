# Plot LZs of non-matching lead-clumped lead loci
library(dplyr)
library(purrr)
library(tidyr)
library(vroom)

source('/Volumes/scratch/merrimanlab/riku/merrimanlab_scratch/metQTL/src/locuszoom/LocusZooms/functions/locus_zoom.R')

genes = read.delim('/Volumes/scratch/merrimanlab/riku/merrimanlab_scratch/metQTL/src/locuszoom/LocusZooms/Gencode_GRCh37_Genes_UniqueList2021.txt', header = T, stringsAsFactors = F)

no_match = read.table('res/urate_clumping/lead_not_match_clump.txt', sep = '\t', header = T, stringsAsFactors = F)
no_match = no_match %>% separate(loci, sep = '_', remove = F, into = c('CHR', 'start', 'end'), convert = T) %>% mutate(CHR = as.numeric(gsub('chr', '', CHR)))

urate = vroom('res/urate_meta/tin_ukbb.urate_meta1.clean.txt')
urate = urate %>% select(CHR, POS, P, SNP)

colnames(urate) = c('CHR', 'BP', 'P', 'SNP')

for (i in 1:nrow(no_match)) {
	snp = no_match$clump_lead[i]
	region = c(no_match$CHR[i], no_match$start[i], no_match$end[i])
	plot_dat = urate %>% filter(CHR == no_match$CHR[i], between(BP, no_match$start[i] - 500000, no_match$end[i] + 500000)) %>% as.data.frame
	out_name = gsub('LOCUS', no_match$loci[i], 'res/urate_clumping/no_clump_lz/LOCUS.png')
	locus.zoom(data = plot_dat,
			   snp = snp,
			   offset_bp = 500000,
			   genes.data = genes,
			   file.name = out_name,
			   ignore.lead = T,
			   rsid.check = F)
}

