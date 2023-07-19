# Make columns with locus boundary information
get_boundary = function(data) {
	res = data %>% mutate(clean_locus = gsub('MB.*', '', locus)) %>% separate(clean_locus, into = c('chr', 'start', 'end'), sep = '_', remove = F, convert = T) %>% mutate(start = start * 1e6, end = end * 1e6)
	res$chr[res$chr == 'chr23'] = 'chrX'
	return(res)
}

# Make initial gene list using locus boundaries
get_overlap_genes = function(data, gencode) {
	data$gene = NA
	for (i in 1:nrow(data)) {
		chr = data$chr[i]
		start_pos = data$start[i]
		end_pos = data$end[i]
		data_gene = gencode %>% filter(Chrom == chr & (between(Start, start_pos, end_pos) | between(End, start_pos, end_pos) | (Start <= start_pos & End >= end_pos)))
		if (nrow(data_gene) > 0) {
			data$gene[i] = paste(data_gene$Gene, collapse = ';')
			next
		}
	}

	# Split the genes into separate rows
	data = data %>% separate_rows(gene, sep = ';')
	data = left_join(data, gencode, by = c('gene' = 'Gene'))
	data = data %>% select(locus, Main_GWAS_loci, ancestry, SNP, CHR, BP, gene:End, ensemblGeneID_short, Coding) %>% distinct
	data$gene_chr = gsub('chr', '', data$Chrom)
	data$gene_chr = gsub('X', '23', data$gene_chr)

	# Since there are some duplicate genes (mainly RNAs) in the gencode file, only
	# keep those that have matching chromosome between the locus and the gencode
	# chromosome
	ind = which(data$CHR == data$gene_chr)
	data = data[ind,]

	colnames(data) = c("locus", "Main_GWAS_loci", "ancestry", "SNP", "snp_chr", "snp_pos", "gene", "Chrom", "gene_start", "gene_end", "ENSG_short", "coding", "gene_chr")

	# Make the initial gene list
	gene_map = data %>% select(locus:gene, gene_chr, gene_start:coding) %>% mutate(within_locus = T, gene_chr = as.integer(gene_chr))

	# Find genes that are completely non-sensical (e.g. Y_RNA, XXbac, SnoU13, etc.)

	# First, list all duplicated genes
	dup_gene = gene_map$gene[which(duplicated(gene_map$gene))]

	# Weird/stupid genes usually have the same gene names but have multiple
	# different ENSG ID:
	bad_gene = gene_map %>% filter(gene %in% dup_gene) %>% group_by(gene, ENSG_short) %>% summarize(count = n()) %>% distinct %>% ungroup %>% select(!count) %>% group_by(gene) %>% summarize(count = n()) %>% filter(count > 1)

	# For each bad gene in the list, determine whether the gene is actually bad or
	# not. Easiest indication is if the gene is "present" in multiple chromosomal
	# location (e.g. in chr1, 2, and 3)
	remove_gene = c()
	for (i in 1:nrow(bad_gene)) {
		tmp = gene_map %>% filter(gene == bad_gene$gene[i])
		n_chr = table(tmp$gene_chr)
		if(length(n_chr) > 1) {
			remove_gene = c(remove_gene, bad_gene$gene[i])
		}
	}

	# Remove the genes that are blatantly wrong
	gene_map = gene_map %>% filter(!(gene %in% remove_gene))

	# Also remove genes that have "XXbac" in their gene name
	gene_map = gene_map %>% filter(!grepl('XXbac', gene))

	# Before moving on, make sure that all the genes are actually within the locus
	# boundary (i.e. remove the remaining weird genes)
	boundary = data %>% mutate(clean_locus = gsub('MB.*', '', locus)) %>% separate(clean_locus, into = c('chr', 'start', 'end'), sep = '_', remove = F, convert = T) %>% mutate(start = start * 1e6, end = end * 1e6) %>% select(locus, snp_chr, ancestry, start:end) %>% distinct

	bad_gene_subset = gene_map %>% filter(gene %in% bad_gene$gene)

	rm_ensg = c()
	for (i in 1:nrow(bad_gene_subset)) {
		tmp = boundary %>% filter(locus == bad_gene_subset$locus[i], ancestry == bad_gene_subset$ancestry[i])
		start_pos = tmp$start
		end_pos = tmp$end
		rm_gene = bad_gene_subset[i,] %>% filter(!(between(gene_start, start_pos, end_pos) | between(gene_end, start_pos, end_pos)))
		if(nrow(rm_gene) > 0) {
			rm_ensg = c(rm_ensg, bad_gene_subset$ENSG_short[i])
		}
	}

	gene_map = gene_map %>% filter(!(ENSG_short %in% rm_ensg))

	return(gene_map)
}

# Clean up eQTL data
clean_eqtl = function(eqtl, data, gencode) {
	# For female eQTL, since ther ewere no eQTL found in female-specific
	# analysis
	if(nrow(eqtl) == 0) {
		return(eqtl)
	}
	eqtl_gene = eqtl %>% mutate(ENSG_short = gsub('\\..*$', '', ENSG)) %>% select(SNP, CHR:BP, gene_name, ENSG_short, locus.type) %>% distinct
	eqtl_gene = left_join(eqtl_gene, data, by = c('SNP'))
	eqtl_gene = eqtl_gene %>% select(locus, ancestry, SNP:locus.type)
	eqtl_gene$type = "eQTL"
	colnames(eqtl_gene)[4:5] = c('snp_chr', 'snp_pos')

	# Some genes from GTEx aren't present in the gencode data, or are present, but
	# with different ENSG ID. Fix some of these:
	alias = read.table('get_alias.txt', sep = '\t', header = F, stringsAsFactors = F)
	colnames(alias) = c('original', 'alias', 'ENSG', 'chrpos')

	alias_sub = alias %>% filter(!is.na(alias))

	for (i in 1:nrow(alias_sub)) {
		ind = which(eqtl_gene$gene_name == alias_sub$original[i])
		eqtl_gene$gene_name[ind] = alias_sub$alias[i]
		eqtl_gene$ENSG_short[ind] = alias_sub$ENSG[i]
	}

	# Add gene chr and start/end positions using gencode data
	eqtl_gene = gencode %>% select(Chrom:End, ensemblGeneID_short, Coding) %>% left_join(eqtl_gene, ., by = c('ENSG_short' = 'ensemblGeneID_short')) %>% select(locus:gene_name, Chrom:End, ENSG_short:Coding)

	# For those genes that I couldn't find an alternative ENSG or alias, use the
	# GRCh37-mapped chrpos taken from ensembl
	alias_sub = alias %>% filter(is.na(alias)) %>% separate(chrpos, into = c('chr', 'start', 'end'), convert = T)

	for (i in 1:nrow(alias_sub)) {
		ind = which(eqtl_gene$gene_name == alias_sub$original[i])
		eqtl_gene$Chrom[ind] = alias_sub$chr[i]
		eqtl_gene$Start[ind] = alias_sub$start[i]
		eqtl_gene$End[ind] = alias_sub$end[i]
	}

	# For the genes that couldn't be found, fill in the 'Coding' column
	proteincoding = c('PCDHA7')
	pseudo = c('RPL23AP97', 'NBPF25P')

	eqtl_gene$Coding[eqtl_gene$gene_name %in% proteincoding] = 'proteincoding'
	eqtl_gene$Coding[eqtl_gene$gene_name %in% pseudo] = 'psuedogene'
	eqtl_gene$Coding[is.na(eqtl_gene$Coding)] = 'lncRNA'

	# Final clean-up
	eqtl_gene$Chrom = gsub('chr', '', eqtl_gene$Chrom)
	eqtl_gene$Chrom = as.integer(gsub('X', '23', eqtl_gene$Chrom))
	colnames(eqtl_gene) = c("locus", "ancestry", "SNP", "snp_chr", "snp_pos", "gene_name", "gene_chr", "gene_start", "gene_end", "ENSG_short", "locus.type", "type", "coding")

	return(eqtl_gene)
}

# Add eQTL gene list to full and male (no female eQTL found):
add_eqtl_gene = function(data, eqtl) {
	if(nrow(eqtl) == 0) {
		data$eqtl = 'NA-NA'
	} else {
		data = full_join(data, eqtl, by = c('locus', 'ancestry', 'SNP', 'snp_chr', 'snp_pos', 'gene' = 'gene_name', 'gene_chr', 'gene_start', 'gene_end', 'ENSG_short', 'coding'))
		data = data %>% unite(col = 'eqtl', locus.type, type, sep = '-')
	}
	# Make an eQTL column and make the within_locus info consistent
	data$eqtl[which(data$eqtl == 'NA-NA')] = F
	data$within_locus[which(is.na(data$within_locus))] = F
	data = as.data.frame(data)
	return(data)
}

