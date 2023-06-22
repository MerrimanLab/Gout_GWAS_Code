# Script to make a table of all the proxies and add missense information
library(dplyr)
library(vroom)
library(tidyr)
library(purrr)

files = list.files('results/missense/', pattern = '.tags.list', full.names = T)

ancestry = gsub('\\.tags.list', '', files)
ancestry = toupper(gsub('.*\\.', '', ancestry))

dat_list = map(files, ~ read.table(.x, sep = '', header = T, stringsAsFactors = F))
dat_list = map(dat_list, ~ .x %>% separate_rows(TAGS))

# There were 30 variants not found in 1KGP data I used for LD, so need to add
# these back in to determine if they are missense:
cand_list = readLines('data/missense/candidate_list.txt')
dat_snp = map(dat_list, ~ .x %>% pull(SNP)) %>% unlist %>% unique
ind = which(!(cand_list %in% dat_snp))
cand_list = data.frame(SNP = cand_list[ind])

# Combine the data, add ancestry info, and remove rs117897057 from the list
dat_list = map2(dat_list, ancestry, ~ full_join(.x, cand_list) %>% mutate(ld_ancestry = .y) %>% filter(SNP != 'rs117897057'))

# Add locus info to the proxy list
loci_list = c('/Volumes/archive/merrimanlab/major_gwas_paper_archive/results_for_paper/loci_list/full_indepSNP_summary_updated_9Dec2022_withBroad.txt', '/Volumes/archive/merrimanlab/major_gwas_paper_archive/results_for_paper/loci_list/male_indepSNP_summary_updated_9Dec2022_withBroad.txt', '/Volumes/archive/merrimanlab/major_gwas_paper_archive/results_for_paper/loci_list/female_indepSNP_summary_updated_9Dec2022_withBroad.txt')
loci = map(loci_list, ~ read.table(.x, sep = '\t', header = T, stringsAsFactors = F))

# Need to adjust/match the male/female-specific loci to the full loci
sync_loci = function(loci){
	ind = grepl('Unique', loci$Main_GWAS_loci)
	loci$Main_GWAS_loci[ind] = loci$locus[ind]
	loci$locus = loci$Main_GWAS_loci
	loci = loci %>% select(!Main_GWAS_loci)
	return(loci)
}

loci[2:3] = map(loci[2:3], ~ sync_loci(.x))
loci = map(loci, ~ .x %>% select(SNP, locus, SEX))
loci = full_join(loci[[1]], loci[[2]], by = c('SNP', 'locus'), suffix = c('', '.male')) %>% full_join(., loci[[3]], by = c('SNP', 'locus'), suffix = c('.full', '.female')) %>% unite(cohort, contains('SEX'), sep = '/', na.rm = T)

dat_list = map(dat_list, ~ left_join(.x, loci))

# For those without a locus, they must have come from the credible set list or COJO list
fm_res = read.table('results/credible_set/combined_results/top_candidate_cred_var.combined.txt', sep = '\t', header = T, stringsAsFactors = F)

dat_list = map(dat_list, ~ fm_res %>% select(locus, SNP) %>% left_join(.x, ., by = c('SNP')) %>% mutate(locus = ifelse(is.na(locus.x), locus.y, locus.x)) %>% select(!c(locus.x, locus.y)))

# Now figure out COJO loci:
cojo_file = list.files('data/cojo_dat', pattern = '.txt', full.names = T)
cojo_list = map2_dfr(cojo_file, c('female', 'male', 'full'), ~ read.table(.x, sep = '\t', header = T, stringsAsFactors = F) %>% mutate(cohort = .y))
cojo_list = cojo_list %>% select(Chr:bp, cohort) %>% pivot_wider(names_from = cohort, values_from = cohort, values_fill = NA) %>% unite(col = 'cohort', full, male, female, sep = '/', na.rm = T)

# Assign locus to COJO variants based on which locus the variant is in:
add_locus = function(data, loci) {
	data$locus = NA
	for (i in 1:nrow(loci)) {
		ind = which(data$Chr == loci$chr[i] & between(data$bp, loci$start[i], loci$end[i]))
		data$locus[ind] = loci$locus[i]
	}
	return(data)
}

tmp_loci = loci %>% mutate(locus_clean = gsub('MB.*', '', locus)) %>% separate(locus_clean, into = c('chr', 'start', 'end'), sep = '_', convert = T) %>% mutate(chr = as.numeric(gsub('chr', '', chr)), locus = gsub('_[[:alpha:]]$', '', locus)) %>% select(locus, chr, start, end) %>% mutate(start = start * 1e6, end = end * 1e6) %>% distinct

cojo_list = add_locus(cojo_list, tmp_loci)

# There are three variants with no locus info, two of which is close enough to
# existing locus to call it correctly (rs13107325 and rs4789698). Remove the
# last variant (rs619057) since this variant is way too far from any locus
# boundaries (at least for now):

cojo_list$locus[cojo_list$SNP == 'rs13107325'] = 'chr4_102.8_103.16MB'
cojo_list$locus[cojo_list$SNP == 'rs4789698'] = 'chr17_80.43_80.53MB'
cojo_list = cojo_list %>% filter(SNP != 'rs619057')

cojo_list = cojo_list %>% select(SNP, locus, cohort)

dat_list = map(dat_list, ~ cojo_list %>% select(SNP, locus) %>% left_join(.x, ., by = c('SNP')) %>% mutate(locus = ifelse(is.na(locus.x), locus.y, locus.x)))

# Remove rs619057:
dat_list = map(dat_list, ~ .x %>% filter(SNP != 'rs619057'))

dat_list = map(dat_list, ~ .x %>% select(locus, SNP, TAGS, ld_ancestry, cohort))
dat = map_dfr(dat_list, ~ .x)

# First figure out all the unique SNP-proxy pairs, then add on which ancestry
# LD it came from
all_proxy = dat %>% filter(!is.na(locus)) %>% select(-c(ld_ancestry, cohort)) %>% distinct
all_proxy = left_join(all_proxy, dat) %>% distinct

all_proxy = all_proxy %>% pivot_wider(names_from = ld_ancestry, values_from = ld_ancestry, values_fill = NA) %>% unite(col = 'ld_ancestry', AFR:TAMA, sep = '/', na.rm = T)
all_proxy$ld_snp = all_proxy$SNP

# Add missense column using the file downloaded from dbSNP:
missense = vroom('data/dbsnp_missense/missense.bed', col_names = c('chr', 'start', 'pos', 'SNP'))
missense$missense = T
missense = missense %>% select(SNP, missense)

all_proxy = left_join(all_proxy, missense, by = c('TAGS' = 'SNP'))
all_proxy$missense[which(is.na(all_proxy$missense))] = F

# Add some info to the variants using a list of all the variants present from
# the gout GWAS from all ancestries
info = vroom('data/missense/all_gwas_variants.txt')
# info = info %>% distinct(SNP, .keep_all = T)

all_proxy = left_join(all_proxy, info, by = c('TAGS' = 'SNP'))

# Check if the initial set of candidate variants are missense
candidate = all_proxy %>% distinct(SNP) %>% left_join(., missense) %>% left_join(., info)
candidate$missense[which(is.na(candidate$missense))] = F
candidate$ld_ancestry = NA
candidate$proxy_of = NA

candidate = candidate %>% select(SNP, CHR:major, missense, proxy_of, ld_ancestry)

# Make a list of all the proxy variants that are missense
proxy_missense = all_proxy %>% filter(!(TAGS %in% candidate$SNP), missense)
proxy_missense = proxy_missense %>% select(ld_snp, TAGS, ld_ancestry, missense:major) %>% distinct
proxy_missense$proxy_of = proxy_missense$ld_snp

proxy_missense = proxy_missense %>% select(TAGS, CHR:major, missense, proxy_of, ld_ancestry)
colnames(proxy_missense)[1] = 'SNP'

# Remove missense variants that are not present in the GWAS data
proxy_missense = proxy_missense %>% filter(!is.na(CHR))

# Bind the two datasets
full_list = rbind(candidate, proxy_missense) %>% distinct

# Add locus and cohort info back
loci_list = c('data/loci_list/full_indepSNP_summary_1jun2022.txt', 'data/loci_list/male_indepSNP_summary_1jun2022.txt', 'data/loci_list/female_indepSNP_summary_1jun2022.txt')
loci = map(loci_list, ~ read.table(.x, sep = '\t', header = T, stringsAsFactors = F))
loci[2:3] = map(loci[2:3], ~ sync_loci(.x))
loci = map(loci, ~ .x %>% select(SNP, locus, SEX:TAMA))
loci = full_join(loci[[1]], loci[[2]], by = c('SNP', 'locus'), suffix = c('', '.male')) %>% full_join(., loci[[3]], by = c('SNP', 'locus'), suffix = c('.full', '.female')) %>% unite(cohort, contains('SEX'), sep = '/', na.rm = T)
loci[loci == ''] = NA
loci[loci == 'LAT (TAMA-LD)'] = 'LAT'
loci$locus_ancestry = apply(loci[, 4:ncol(loci)], 1, function(x) paste(unique(x[!is.na(x)]), collapse = '/'))
loci = loci %>% select(SNP:cohort, locus_ancestry) %>% mutate(locus = gsub('_[[:alpha:]]$', '', locus))
loci$variant_from = 'lead'

tmp_fm = fm_res %>%
	select(locus, SNP, cohort) %>%
	distinct %>%
	pivot_wider(names_from = cohort, values_from = cohort, values_fill = NA) %>%
	unite(col = 'cohort', full, male, female, sep = '/', na.rm = T) %>%
	distinct %>%
	mutate(variant_from = 'finemapping')

# Add cohort/loci info from finemapping
full_list = left_join(full_list, tmp_fm, by = 'SNP') %>% left_join(., tmp_fm, by = c('proxy_of' = 'SNP'))
full_list = full_list %>% mutate(locus = ifelse(is.na(locus.y), locus.x, locus.y), cohort = ifelse(is.na(cohort.y), cohort.x, cohort.y)) %>% unite(variant_from, contains('variant_from'), na.rm = T, sep = '/') %>% select(locus, cohort, SNP:ld_ancestry, variant_from)
full_list$variant_from[full_list$variant_from == ''] = NA

# Add info from loci list (information from loci list will be prioritised over
# the finemapping info in terms of locus/cohort):
full_list = left_join(full_list, loci, by = c('SNP')) %>% left_join(., loci, by = c('proxy_of' = 'SNP'))
full_list = full_list %>% mutate( locus.y = ifelse(is.na(locus), locus.y, locus), cohort.y = ifelse(is.na(cohort), cohort.y, cohort), locus_ancestry = ifelse(is.na(locus_ancestry.x), locus_ancestry.y, locus_ancestry.x), variant_from.y = ifelse(is.na(variant_from), variant_from.y, variant_from)) %>% select(locus.x:variant_from.y, locus_ancestry)

full_list = full_list %>% mutate(locus = ifelse(is.na(locus.y), locus.x, locus.y), cohort = ifelse(is.na(cohort.y), cohort.x, cohort.y)) %>% unite(variant_from, contains('variant_from'), na.rm = T, sep = '/') %>% select(locus, cohort, locus_ancestry, SNP:ld_ancestry, variant_from)
full_list$variant_from[full_list$variant_from == ''] = NA

# Add COJO info:
full_list = cojo_list %>% mutate(variant_from = 'COJO') %>% left_join(full_list, ., by = 'SNP')
full_list = cojo_list %>% mutate(variant_from = 'COJO') %>% left_join(full_list, ., by = c('proxy_of' = 'SNP'))
full_list = full_list %>% mutate( locus.y = ifelse(is.na(locus), locus.y, locus), cohort.y = ifelse(is.na(cohort), cohort.y, cohort), variant_from.y = ifelse(is.na(variant_from), variant_from.y, variant_from)) %>% select(locus.x:variant_from.y)

full_list = full_list %>% mutate(locus = ifelse(is.na(locus.x), locus.y, locus.x), cohort = ifelse(is.na(cohort.x), cohort.y, cohort.x)) %>% unite(variant_from, contains('variant_from'), na.rm = T, sep = '/') %>% select(locus, locus_ancestry, cohort, SNP:ld_ancestry, variant_from)

# Final clean up:
# - add ancestry from finemapping
tmp_fm2 = fm_res %>% select(SNP, contains('pip'))
tmp_fm2$locus_ancestry = ifelse(is.na(tmp_fm2$pip.tama), 'EUR', ifelse(is.na(tmp_fm2$pip.paintor) & is.na(tmp_fm2$pip.finemap), 'TAMA', 'EUR/TAMA'))
tmp_fm2 = tmp_fm2 %>% select(SNP, locus_ancestry)

full_list = left_join(full_list, tmp_fm2, by = 'SNP')
full_list = left_join(full_list, tmp_fm2, by = c('proxy_of' = 'SNP'))

full_list = full_list %>% mutate(locus_ancestry.y = ifelse(is.na(locus_ancestry), locus_ancestry.y, locus_ancestry)) %>% select(locus:locus_ancestry.y)
full_list = full_list %>% mutate(locus_ancestry = ifelse(is.na(locus_ancestry.x), locus_ancestry.y, locus_ancestry.x)) %>% select(locus, locus_ancestry, cohort:variant_from)

# - any NA ancestry that are from COJO will be 'EUR'
full_list = cojo_list %>% mutate(locus_ancestry = 'EUR') %>% distinct(SNP, locus_ancestry) %>% left_join(full_list, ., by = 'SNP')
full_list = cojo_list %>% mutate(locus_ancestry = 'EUR') %>% distinct(SNP, locus_ancestry) %>% left_join(full_list, ., by = c('proxy_of' = 'SNP'))

full_list = full_list %>% mutate(locus_ancestry.y = ifelse(is.na(locus_ancestry), locus_ancestry.y, locus_ancestry)) %>% select(locus:locus_ancestry.y)
full_list = full_list %>% mutate(locus_ancestry = ifelse(is.na(locus_ancestry.x), locus_ancestry.y, locus_ancestry.x)) %>% select(locus, locus_ancestry, cohort:variant_from)

# Save the full list and the list of missense variants:
missense_snp = full_list %>% filter(missense) %>% pull(SNP) %>% unique
writeLines(missense_snp, 'results/missense/missense_variants.txt')
write.table(full_list, 'results/missense/candidate_variants.missense_info.txt', sep = '\t', col.names = T, row.names = F, quote = F)

