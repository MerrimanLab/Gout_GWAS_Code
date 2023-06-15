# Script to generate the finemapping stats for the paper
library(dplyr)
library(purrr)
library(tidyr)
library(ggplot2)
library(Manu)

################################################################################
# Load the list of variants with PIP >= 0.5 - use this to pull out the variants
# in the credible set from the raw data
cohort = c('full', 'male', 'female')
comb_file = map(cohort, ~ gsub('sex', .x, 'results/credible_set/combined_results/sex/combined_finemap_result.sex.txt'))
dat_list = map2(comb_file, cohort, ~ read.table(.x, sep = '\t', header = T, stringsAsFactors = F) %>% mutate(cohort = .y))

# Load all the raw finemapping data
paintor_file = map(cohort, ~ gsub('sex', .x, 'results/credible_set/PAINTOR/output/sex/credible_pip_list.sex.txt'))
multi_file   = map(cohort, ~ gsub('sex', .x, 'results/credible_set/PAINTOR/output_multi/sex/credible_pip_list.sex.txt'))
finemap_file = map(cohort, ~ gsub('sex', .x, 'results/credible_set/FINEMAP/sss_default/sex/credible_table.sex.txt'))
tama_file    = map(cohort, ~ gsub('sex', .x, 'results/credible_set/TAMA/sex/tama_sex_abf.credible_list.txt'))

paintor       = map2(paintor_file , cohort, ~ read.table(.x, sep = '\t', header = T, stringsAsFactors = F) %>% mutate(cohort = .y))
paintor_multi = map2(multi_file   , cohort, ~ read.table(.x, sep = '\t', header = T, stringsAsFactors = F) %>% mutate(cohort = .y))
finemap       = map2(finemap_file , cohort, ~ read.table(.x, sep = '\t', header = T, stringsAsFactors = F) %>% mutate(cohort = .y))
tama          = map2(tama_file    , cohort, ~ read.table(.x, sep = '\t', header = T, stringsAsFactors = F) %>% mutate(cohort = .y))

################################################################################
# For PAINTOR results, replace some loci with the multiple causal variant
# results:
paintor = map2(paintor, paintor_multi, ~ .x %>% filter(!(locus.SNP %in% .y$locus.SNP)))
paintor = map2(paintor, paintor_multi, ~ rbind(.x, .y))

# Merge locus info to PAINTOR data using FINEMAP data
tmp_loci = map(finemap, ~ .x %>% select(locus, SNP.lead) %>% distinct)
paintor = map2(paintor, tmp_loci, ~ left_join(.x, .y, by = c('locus.SNP' = 'SNP.lead')) %>% rename(pip = Posterior_Prob))

colnames(paintor)[10] = 'pip'

# Determine whether the credible set has variant with PIP >= 0.5
check_pip = function(data) {
	tmp_dat = data %>%
		group_by(locus) %>%
		arrange(desc(pip)) %>%
		slice(1) %>%
		ungroup %>%
		mutate(good = pip >= 0.5) %>%
		select(locus, good)
	return(tmp_dat)
}

tmp_pt = map(paintor, check_pip)
paintor = map2(paintor, tmp_pt, ~ left_join(.x, .y))

################################################################################
# Add locus info to the TAMA results:
loci_file = map(cohort, ~ gsub('sex', .x, '/Volumes/archive/merrimanlab/major_gwas_paper_archive/results_for_paper/loci_list/sex_indepSNP_summary_updated_9Dec2022_withBroad.txt'))
tama_loci = map(loci_file, ~ read.table(.x, sep = '\t', header = T, stringsAsFactors = F) %>% filter(TAMA == 'TAMA') %>% rename(cohort = SEX))

clean_tama_loci = function(loci) {
	tmp_loci = loci
	tmp_loci$locus_orig = tmp_loci$locus
	tmp_loci$locus = gsub('[[:alpha:]]', '', tmp_loci$locus)
	tmp_loci$locus = gsub('_$', '', tmp_loci$locus)
	tmp_loci = tmp_loci %>% separate(locus, sep = '_', into = c('CHR', 'start', 'end'), remove = F, convert = T)
	tmp_loci$locus = tmp_loci$locus_orig
	tmp_loci = tmp_loci %>% select(locus:end, SNP, cohort) %>% mutate(start = start * 1000000, end = end * 1000000)
	return(tmp_loci)
}

tama_loci = map(tama_loci, clean_tama_loci)

# Pull out TAMA finemapping results using the loci region:
pull_tama_regions = function(data, loci) {
	tmp_dat = data %>% mutate(locus = NA, locus.SNP = NA) %>% select(locus, locus.SNP, SNP, CHR, POS, post_prob, cohort)
	for (i in 1:nrow(loci)) {
		chr = loci$CHR[i]
		start = loci$start[i]
		end = loci$end[i]
		ind = which(tmp_dat$CHR == chr & between(tmp_dat$POS, start, end))
		tmp_dat$locus[ind] = loci$locus[i]
		tmp_dat$locus.SNP[ind] = loci$SNP[i]
	}
	colnames(tmp_dat)[6] = 'pip'
	return(tmp_dat)
}

tama = map2(tama, tama_loci, ~ pull_tama_regions(.x, .y))

# Determine whether the credible set has variant with PIP >= 0.5
tmp_bf = map(tama, check_pip)
tama = map2(tama, tmp_bf, ~ left_join(.x, .y))

################################################################################
# Clean up FINEMAP data so that it only includes variants in a credible set
# with at least one variant with PIP >= 0.5
finemap = map(finemap, ~ .x %>% select(locus, SNP.lead, credible_set:pip, cohort))

# Determine whether the credible set has variant with PIP >= 0.5.
# Note that I can't re-use the other function since FINEMAP results require
# different grouping
for (i in 1:length(finemap)) {
	tmp_fm = finemap[[i]] %>%
	group_by(locus, credible_set) %>%
	arrange(desc(pip)) %>%
	slice(1) %>%
	ungroup %>%
	mutate(good = pip >= 0.5) %>%
	select(locus:total_snp_in_cs, good)
	# Merge this information:
	finemap[[i]] = left_join(finemap[[i]], tmp_fm)
}

################################################################################
# Need to sync up the loci so full/female/male loci matches up
# - merge "main gwas loci" to male/female
# - add male/female/full column
# - rbind and group by locus

loci_list = map(loci_file[2:3], ~ read.table(.x, sep = '\t', header = T, stringsAsFactors = F) %>% select(locus, Main_GWAS_loci) %>% distinct)

# Edit the sex_specific locus so that it matches the loci listed in the full
# cohort loci list
sync_loci = function(data, loci){
	data = left_join(data, loci)
	ind = grepl('Unique', data$Main_GWAS_loci)
	data$Main_GWAS_loci[ind] = data$locus[ind]
	data$locus = data$Main_GWAS_loci
	data = data %>% select(!Main_GWAS_loci)
	return(data)
}

finemap[2:3] = map2(finemap[2:3], loci_list, ~ sync_loci(.x, .y))
tama[2:3] = map2(tama[2:3], loci_list, ~ sync_loci(.x, .y))
paintor[2:3] = map2(paintor[2:3], loci_list, ~ sync_loci(.x, .y))

finemap = map_dfr(finemap, ~ .x %>% mutate(locus = gsub('_[[:alpha:]]$', '', locus))) %>% mutate(program = 'FINEMAP') %>% rename(locus.SNP = SNP.lead)
tama = map_dfr(tama, ~ .x %>% mutate(locus = gsub('_[[:alpha:]]$', '', locus))) %>% mutate(program = 'TAMA_BF')
paintor = map_dfr(paintor, ~ .x %>% mutate(locus = gsub('_[[:alpha:]]$', '', locus))) %>% mutate(program = 'PAINTOR')

# Count the total unique loci that were finemapped from EUR, TAMA, and both:
eur_uniq_loci = unique(c(finemap$locus, paintor$locus))
tama_uniq_loci = unique(tama$locus)
total_uniq_loci = unique(c(eur_uniq_loci, tama_uniq_loci))

length(eur_uniq_loci)
length(tama_uniq_loci)
length(total_uniq_loci)

################################################################################
# Merge all the variants from all the finemapping programs
res = paintor %>% select(locus, locus.SNP:SNP, pip, good, cohort)
res = finemap %>% select(!program) %>% full_join(res, ., by = c('locus', 'locus.SNP', 'SNP' = 'SNP.credible', 'cohort'), suffix = c('', '.finemap'))
res = tama %>% select(locus, locus.SNP, SNP, pip, good, cohort) %>% full_join(res, ., by = c('locus', 'locus.SNP', 'SNP', 'cohort'), suffix = c('.paintor', '.tama'))

# Pull out credible variants from loci with at least one variant with PIP >=
# 0.5 in any one of the three finemapping programs:
res$good_any = pmap_lgl(list(res$good.paintor, res$good.tama, res$good.finemap), ~ any(..1, ..2, ..3, na.rm = T))

res = res %>% filter(good_any) %>% select(!contains('good')) %>% arrange(locus)

# Count the number of loci with "good" credible set:
length(unique(res$locus))

################################################################################
# Count up how many of the loci are suspicious or not.

# Sync the locus naming of SLALOM results:
dat_list[2:3] = map2(dat_list[2:3], loci_list, ~ sync_loci(.x, .y))

# A locus is deemed suspicious if either EUR or TAMA SLALOM result flag the
# locus as suspicious. This is probably more conservative than it should be,
# but I think it's okay, considering there will be a LOT of false positives
# from fine-mapping.
suspicious = map(dat_list, ~ .x %>% filter(suspicious_slalom.eur | suspicious_slalom.tama) %>% pull(locus) %>% unique)
suspicious = map(suspicious, ~ gsub('_[[:alpha:]]$', '', .x) %>% unique)

tmp = map(dat_list, ~ .x %>% pull(locus) %>% gsub('_[[:alpha:]]$', '', .) %>% unique)

length(unique(unlist(suspicious))) # 101 suspicous loci in total

# How many of the 285 loci are suspicious?
length(which(unique(res$locus) %in% unique(unlist(suspicious)))) # 101

# For each cohort, remove the suspicious loci
not_susp = map2_dfr(cohort, suspicious, ~ res %>% filter(cohort == .x) %>% filter(!(locus %in% .y)))

length(unique(not_susp$locus)) # 215 unique loci remaining from various cohorts

################################################################################
# Pull out credible sets with <= 5 variants in >1 fine-mapping methods
nonsusp_loci = map(cohort, ~ not_susp %>% filter(cohort == .x) %>% pull(locus) %>% unique)

finemap_nonsusp = map2_dfr(cohort, nonsusp_loci, ~ finemap %>% filter(cohort == .x) %>% filter(locus %in% .y & good))
tama_nonsusp = map2_dfr(cohort, nonsusp_loci, ~ tama %>% filter(cohort == .x) %>% filter(locus %in% .y & good))
paintor_nonsusp = map2_dfr(cohort, nonsusp_loci, ~ paintor %>% filter(cohort == .x) %>% filter(locus %in% .y & good))

# Pull out variants included in credible sets with <= 5 variants from each
# approach
tama_cv = tama_nonsusp %>% group_by(locus, cohort) %>% summarize(count = n()) %>% ungroup %>% left_join(tama_nonsusp, .) %>% filter(`count` <= 5 & pip >= 0.5)
paintor_cv = paintor_nonsusp %>% group_by(locus, cohort) %>% summarize(count = n()) %>% ungroup %>% left_join(paintor_nonsusp, .) %>% filter(`count` <= 5 & pip >= 0.5)
finemap_cv = finemap_nonsusp %>% group_by(locus, cohort) %>% summarize(count = n()) %>% ungroup %>% left_join(finemap_nonsusp, .) %>% filter(`count` <= 5 & pip >= 0.5)

# List out and count the unique variants found by all three approaches
#
# Variant can be listed more than twice if there are two signals at the same
# locus or it appears in different cohorts
all_candidate = c(unique(tama_cv$SNP), unique(paintor_cv$SNP), unique(finemap_cv$SNP.credible))
top_candidate = unique(all_candidate)

length(top_candidate)
table(table(all_candidate))

# Pull out the per-cohort candidate variants from the full table:
candidate = finemap_cv %>% rename(SNP = SNP.credible) %>% list(., tama_cv, paintor_cv) %>% map_dfr(., ~ .x %>% select(locus.SNP, SNP, cohort)) %>% distinct %>% mutate(candidate = T)

loci = map(loci_file, ~ read.table(.x, sep = '\t', header = T, stringsAsFactors = F))
loci[2:3] = map2(loci[2:3], loci[2:3], ~ sync_loci(.x, .y))
loci = map_dfr(loci, ~ .x %>% select(locus: SNP)) %>% distinct

top_res = left_join(not_susp, candidate)
top_res$candidate[is.na(top_res$candidate)] = F
top_res = top_res %>% filter(candidate) %>% select(!candidate) %>% distinct %>% select(locus:SNP, cohort, contains('pip'))
top_res = left_join(top_res, loci, by = c('locus', 'locus.SNP' = 'SNP')) %>% select(locus, CHR:BP, locus.SNP:pip.tama) %>% arrange(CHR, BP, cohort) %>% select(-c(CHR, BP))

write.table(top_res,  'results/credible_set/combined_results/top_candidate_cred_var.combined.txt', sep = '\t', col.names = T, row.names = F, quote = F)

