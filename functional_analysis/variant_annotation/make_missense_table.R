# Script to merge the missense variant annotation, LD information, and gout
# risk allele
#
# Run this script for full/male/female and then combine the three using src/missense/combine_missense_table.R
library(dplyr)
library(tidyr)
library(stringr)
library(purrr)

args = commandArgs(trailingOnly = T)

cohort = args[1]

# Load missense variant annotation:
annot = read.delim('results/missense/vep_output.txt', sep = '\t', header = T, stringsAsFactors = F, comment.char = '')
colnames(annot)[1] = 'SNP'

annot = annot %>% filter(BIOTYPE == 'protein_coding', EXON != '-', grepl('missense', Consequence)) %>% distinct

annot_small = annot %>% select(SNP:Gene, BIOTYPE, Protein_position:Codons, CADD_PHRED, STRAND) %>% separate(Location, into = c('CHR', 'POS', 'POS1'), sep = '[:-]') %>% select(!POS1) %>% distinct

# Figure out the proxy variants of the missense variants (if any)
proxy = read.table('results/missense/missense_gene.txt', sep = '\t', header = T, stringsAsFactors = F)
proxy = proxy %>% select(!c(gene, missense))

# Assign alleleic value to check whether the GWAS variant matches the missense
# variant later on:
minor_val = unlist(lapply(proxy$minor, function(x) switch(x, A = 1, C = 2, G = 3, T = 5)))
major_val = unlist(lapply(proxy$major, function(x) switch(x, A = 1, C = 2, G = 3, T = 5)))
proxy$av = minor_val + major_val

# Join the missense annotations and remove variants that don't have annotations
# (likely due to only synonymous mutations present from VEP):
# (NOTE: 8 variants get kicked out)
res = left_join(proxy, annot_small, by = 'SNP', suffix = c('', '.annot'), keep = T)
length(unique(res$SNP))
res = res %>% filter(!is.na(SNP.annot))
length(unique(res$SNP))

# Split the codons and take out the alleles
res = res %>% separate(Codons, sep = '/', into = c('codon1', 'codon2'), remove = F) %>% mutate(miss_a1 = str_remove_all(codon1, '[[:lower:]]'), miss_a2 = str_remove_all(codon2, '[[:lower:]]'))

# Assign allelic value to the codon changing alleles
minor_val = unlist(lapply(res$miss_a1, function(x) switch(x, A = 1, C = 2, G = 3, T = 5)))
major_val = unlist(lapply(res$miss_a2, function(x) switch(x, A = 1, C = 2, G = 3, T = 5)))
res$av1 = minor_val + major_val

# Flip the codon changing alleles to the other strand and assign allelic value:
res$opp_a1 = unlist(lapply(res$miss_a1, function(x) switch(x, A = 'T', C = 'G', G = 'C', T = 'A')))
res$opp_a2 = unlist(lapply(res$miss_a2, function(x) switch(x, A = 'T', C = 'G', G = 'C', T = 'A')))

minor_val = unlist(lapply(res$opp_a1, function(x) switch(x, A = 1, C = 2, G = 3, T = 5)))
major_val = unlist(lapply(res$opp_a2, function(x) switch(x, A = 1, C = 2, G = 3, T = 5)))
res$av2 = minor_val + major_val

# Match either allelic values with the GWAS allelic value.
# If either doesn't match, then it means the missense variant is multi-allelic
# and the annotation must be for a different allele set.
# (220 unique missense variants leftover)
res = res %>% filter((av == av1 | av == av2))
length(unique(res$SNP))

# Clean up a little
res = res %>% select(SNP:ld_ancestry, Allele:Codons, CADD_PHRED:STRAND)

# Add LD information
ld_list = list.files('results/missense/', pattern = '.clean.ld', full.names = T, include.dirs = T)

ancestry = gsub('.clean.ld', '', ld_list)
ancestry = toupper(gsub('.*\\.', '', ancestry))

ld_info = map(ld_list, ~ read.delim(.x, sep = '', header = T, stringsAsFactors = F) %>% filter(SNP_A != SNP_B, SNP_A %in% res$proxy_of))
ld_info = map2_dfr(ld_info, ancestry, ~ .x %>% mutate(ancestry = .y))

ld_info = ld_info %>% spread(key = ancestry, value = R2, sep = '.')
colnames(ld_info) = gsub('ancestry', 'R2', colnames(ld_info))

res = ld_info %>% select(SNP_A, SNP_B, contains('R2')) %>% left_join(res, ., by = c('SNP' = 'SNP_B', 'proxy_of' = 'SNP_A'))

# Load missense variant GWAS summary stats
file_list = paste('results/missense/missense_variants.', ancestry, '_info.', cohort, '.txt', sep = '')

sum_stats = map(file_list, ~ read.table(.x, header = T, sep = '\t', stringsAsFactors = F))

# Focus on certain columns and change the column names just for the TAMA stats:
sum_stats[[5]] = sum_stats[[5]] %>% select(SNP:OTH, weighted_effect)
sum_stats[[5]]$MAF = NA
colnames(sum_stats[[5]]) = c('SNP', 'CHR', 'BP', 'minor', 'major', 'effect', 'MAF')

tmp = map2_dfr(sum_stats, ancestry, ~ .x %>% select(SNP, CHR, BP, minor, major, effect, MAF) %>% mutate(ancestry = .y, minor_is_risk = effect > 0))

sum_stats = tmp %>% group_by(SNP) %>% summarise(count_minor_risk = sum(minor_is_risk)) %>% ungroup %>% left_join(tmp %>% select(!minor_is_risk) %>% pivot_wider(names_from = ancestry, values_from = c('effect', 'MAF'), names_sep = '.'), .) %>% as.data.frame
sum_stats$likely_risk_allele = ifelse(sum_stats$count_minor_risk > 2, sum_stats$minor, sum_stats$major)
sum_stats = sum_stats %>% select(!count_minor_risk)

res = res %>% left_join(., sum_stats, by = c('SNP', 'CHR', 'BP', 'minor', 'major')) %>% select(SNP:major, Allele, likely_risk_allele, contains('effect'), contains('MAF'), proxy_of:ld_ancestry, contains('R2'), Consequence:STRAND)
colnames(res)[6] = 'protein_changing_allele'

# Add posterior probabilities (for those variants that have one)
finemap_files = list.files(gsub('ZZZ', cohort, 'results/credible_set/FINEMAP/sss_default/ZZZ/'), pattern = 'snp', full.names = T, include.dirs = T)
finemap = map_dfr(finemap_files, ~ read.table(.x, sep = ' ', header = T, stringsAsFactors = F)) %>% select(!index) %>% distinct
finemap = finemap %>% select(rsid:position, prob)

tama = read.table(gsub('ZZZ', cohort, 'results/credible_set/TAMA/ZZZ/tama_ZZZ_abf.full_pip.txt'), sep = '\t', header = T, stringsAsFactors = F)
tama = tama %>% select(SNP:POS, post_prob)

paintor_files = list.files(gsub('ZZZ', cohort, 'results/credible_set/PAINTOR/output/ZZZ/'), pattern = 'rs.*results', full.names = T, include.dirs = T)
paintor = map_dfr(paintor_files, ~ read.table(.x, sep = ' ', header = T, stringsAsFactors = F)) %>% distinct
paintor = paintor %>% select(CHR, BP, Posterior_Prob)

# Load in the multi signal PAINTOR output just in case
paintor_multi_files = list.files(gsub('ZZZ', cohort, 'results/credible_set/PAINTOR/output_multi/ZZZ/'), pattern = 'rs.*results', full.names = T, include.dirs = T)
paintor_multi = map_dfr(paintor_multi_files, ~ read.table(.x, sep = ' ', header = T, stringsAsFactors = F)) %>% distinct
paintor_multi = paintor_multi %>% select(CHR, BP, Posterior_Prob)

# Check (very crudely) if any of the missense variants are in the multi-signal loci:
# (None are present, so ignore them)
length(which(res$POS %in% paintor_multi$BP))

# Merge FINEMAP, PAINTOR, and TAMA posterior probabilities:
post = left_join(tama, finemap, by = c('SNP' = 'rsid', 'CHR' = 'chromosome', 'POS' = 'position'))
post = left_join(post, paintor, by = c('CHR', 'POS' = 'BP'))

colnames(post)[4:6] = paste('post_prob', c('TAMA', 'FINEMAP', 'PAINTOR'), sep = '.')

post_subset = post %>% filter(SNP %in% res$SNP) %>% distinct

# There are some duplicate rows due to slightly different posterior probabities
# that came from fine-mapping different signal at the same loci. Most are very
# low PP, so just remove the duplicated rows:
post_subset = post_subset[!duplicated(post_subset$SNP),]

# Add the posterior probabilities to the annotation file:
res = left_join(res, post_subset, by = c('SNP', 'CHR', 'BP' = 'POS'))

res$max_post_prob = pmax(res$post_prob.TAMA, res$post_prob.FINEMAP, res$post_prob.PAINTOR)

# Combine the protein positions for the same missense variant at different
# position (due to variant in different transcripts):
res = res %>% group_by(SNP, proxy_of, Codons) %>% mutate(Protein_position = paste(Protein_position, collapse = '/')) %>% ungroup %>% distinct %>% as.data.frame

# Clean up result:
res = res %>% select(SNP:major, SYMBOL, protein_changing_allele, Protein_position:Codons, STRAND, likely_risk_allele:R2.TAMA, post_prob.TAMA:max_post_prob, CADD_PHRED) %>% filter(CADD_PHRED != '-')

# Save the results:
out_file = gsub('ZZZ', cohort, 'results/missense/missense_variants_information.ZZZ.txt')
write.table(res, out_file, sep = '\t', col.names = T, row.names = F, quote = F)

