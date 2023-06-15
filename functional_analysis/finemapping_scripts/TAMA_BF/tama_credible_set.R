# Script to generate credible sets based on the final list of loci
library(vroom)
library(dplyr)
library(tidyr)
library(purrr)
library(Brobdingnag)

args = commandArgs(trailingOnly = T)

dat = vroom(args[1])

sex = args[2]

loci = read.table(args[3], sep = '\t', header = T, stringsAsFactors = F)
loci = loci %>% filter(TAMA == 'TAMA')

# Make region based on the loci boundaries
loci$locus = gsub('MB.*', '', loci$locus)
loci$locus = gsub('chr', '', loci$locus)

loci = loci %>% separate(locus, into = c('chr', 'start', 'end'), sep = '_', remove = F)
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

# Function to calculate posterior probability of a region
calc_post_prob = function(data, loci) {
	result = list()
	for (i in 1:nrow(loci)) {
		locus = loci$locus[i]
		locus_snp = loci$SNP[i]
		chr = loci$CHR[i]
		start = loci$start[i]
		end = loci$end[i]
		tmp = data %>% filter(CHR == chr, between(POS, start, end)) %>% arrange(desc(logBF))
		# Calculate total regional BF and posterior probability for each SNP
		total_bf = sum(10 ^ as.brob(tmp$logBF))
		tmp$post_prob = as.numeric((10 ^ as.brob(tmp$logBF)) / total_bf)
		tmp$locus = locus
		tmp$locus.SNP = locus_snp
		result[[i]] = tmp
	}
	return(result)
}

# Generate credible set, given a data frame with posterior probability
get_credible = function(data) {
	data = data %>% arrange(desc(post_prob))
	credible = c()
	credible_ppa = 0
	i = 1
	while(credible_ppa < 0.99) {
		credible = c(credible, data$SNP[i])
		credible_ppa = credible_ppa + data$post_prob[i]
		i = i + 1
	}
	credible_list = data[which(data$SNP %in% credible), ]
	return(credible_list)
}

# Function to summarise credible data
summarise_credible = function(credible_list, loci) {
	results = loci %>% select(locus:SNP)
	results$credible_start = map_dbl(credible_list, ~ min(.x$POS))
	results$credible_end = map_dbl(credible_list, ~ max(.x$POS))
	results$credible_range = results$credible_end - results$credible_start
	results$credible_nsnp = map_dbl(credible_list, ~ nrow(.x))
	results$pip_snp_0.1 = map_chr(credible_list, ~ paste(.x$SNP[.x$post_prob >= 0.1], collapse = ','))
	results$pip_snp_0.5 = map_chr(credible_list, ~ paste(.x$SNP[.x$post_prob >= 0.5], collapse = ','))
	return(results)
}

# Make a list of finemapping regions with posterior probability
loci_list = calc_post_prob(dat, loci)

# Save this list:
full_list = map_dfr(loci_list, ~ .x) %>% distinct %>% arrange(CHR, POS)

output = gsub('sex', sex, 'results/credible_set/TAMA/sex/tama_sex_abf.full_pip.txt')
write.table(full_list, output, sep = '\t', col.names = T, row.names = F, quote = F)

# Make a list of 99% credible sets
credible_list = map(loci_list, get_credible)

# Make a list of data frame with credible set information
credible_summary = summarise_credible(credible_list, loci)

# Combine each of the 99% credible sets into separate data frames
credible_data = map_dfr(credible_list, ~ .x %>% select(locus, locus.SNP, SNP:OTH, logBF, post_prob, contains('effect')))

# Also make a list of credible variants with posterior probability >= 0.1 and 0.5
credible_data_0.1 = credible_data %>% filter(post_prob >= 0.1)
credible_data_0.5 = credible_data %>% filter(post_prob >= 0.5)

# Save output
write.table(credible_summary, gsub('full_pip', 'credible_summary', output), sep = '\t', col.names = T, row.names = F, quote = F)
write.table(credible_data, gsub('full_pip', 'credible_list', output), sep = '\t', col.names = T, row.names = F, quote = F)
write.table(credible_data_0.1, gsub('full_pip', 'credible_var.pip0.1', output), sep = '\t', col.names = T, row.names = F, quote = F)
write.table(credible_data_0.5, gsub('full_pip', 'credible_var.pip0.5', output), sep = '\t', col.names = T, row.names = F, quote = F)
