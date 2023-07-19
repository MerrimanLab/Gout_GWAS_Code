library(dplyr)
library(tidyr)
library(purrr)
library(vroom)

source('gene_prioritisation_functions.R')

# Load GTEx eQTL:
eqtl = read.table('/Volumes/archive/merrimanlab/major_gwas_paper_archive/results_for_paper/eQTL/eqtl_results.effects_aligned.final.pp0.8.txt', sep = '\t', header = T, stringsAsFactors = F)

# Load in Hongbo/Suzstak kidney eQTL results:
tub = read.table('/Volumes/scratch/merrimanlab/major_gwas_paper_scratch/rikutakei/revision/dat/kidney_eqtl/Tub.eQTL.coloc.1M.H4.0.5.txt', sep = '\t', header = T, stringsAsFactors = F)
glom = read.table('/Volumes/scratch/merrimanlab/major_gwas_paper_scratch/rikutakei/revision/dat/kidney_eqtl/Glom.eQTL.coloc.1M.H4.0.5.txt', sep = '\t', header = T, stringsAsFactors = F)

tub = tub %>% filter(PP.H4.abf >= 0.8) %>% mutate(tissue = 'Kidney_tubule', source = 'Susztak')
glom = glom %>% filter(PP.H4.abf >= 0.8) %>% mutate(tissue = 'Kidney_glomerular', source = 'Susztak')

colnames(tub) = gsub('Tub_', '', colnames(tub))
colnames(glom) = gsub('Glom_', '', colnames(glom))

kidney = rbind(tub, glom) %>% mutate(locus.type = 'cis', MAF.kidney = NA, cohort = 'full')

# Combine kidney eQTL:
kid_sub = kidney %>% mutate(CHR = as.integer(gsub(':.*', '', SNP))) %>% select(IndexRSID:RSID, CHR, POS:GWAS_Major, tissue, GeneSymbol, Gene, locus.type, GWAS_MAF:GWAS_PVAL, MAF.kidney, eQTL_BETA:eQTL_PVAL, PP.H4.abf, cohort, source)

colnames(kid_sub) = c("SNP", "proxy.SNP", "CHR", "BP", "minor", "major", "tissue", "gene_name", "ENSG", "locus.type", "MAF.gwas", "effect.gwas", "SE.gwas", "P.gwas", "MAF.tissue", "effect.tissue", "SE.tissue", "P.tissue", "PP.H4.abf", "cohort", "source")

# Load meQTL
meqtl = read.table('/Volumes/archive/merrimanlab/major_gwas_paper_archive/results_for_paper/meQTL/coloc_results/meqtl_coloc.pp0.8_n100.all.txt', sep = '\t', header = T, stringsAsFactors = F)

gencode = read.table('/Volumes/scratch/merrimanlab/major_gwas_paper_scratch/rikutakei/GWAS_Functional/data/gene_prioritisation/Gencode_GRCh37_Genes_UniqueList2021.txt', sep = '\t', header = T, stringsAsFactors = F)
gencode$ensemblGeneID = gsub('_[0-9]*$', '', gencode$ensemblGeneID)
gencode$ensemblGeneID_short = gsub('\\..*$', '', gencode$ensemblGeneID)

# Use loci from all the ancestries
loci_file = map_chr(c('full', 'male', 'female'), ~ gsub('sex', .x, '/Volumes/archive/merrimanlab/major_gwas_paper_archive/results_for_paper/loci_list/sex_indepSNP_summary_updated_9Dec2022_withBroad.txt'))
loci = map(loci_file, ~ read.table(.x, sep = '\t', header = T, stringsAsFactors = F))
loci = map(loci, ~ .x %>% select(locus:SNP, AFR:TAMA, SEX))
names(loci) = c('full', 'male', 'female')

loci = map(loci, ~ .x %>% select(locus:TAMA) %>% pivot_longer(AFR:TAMA, names_to = 'dummy', values_to = 'ancestry') %>% filter(ancestry != '') %>% select(locus:SNP, ancestry))
loci[[1]] = loci[[1]] %>% mutate(Main_GWAS_loci = locus) %>% select(locus, Main_GWAS_loci, CHR: ancestry)

# There are some LAT loci that is labelled "LAT (TAMA-LD)" - convert these into
# just "LAT"
loci[[3]]$ancestry[grepl('LAT', loci[[3]]$ancestry)] = 'LAT'

################################################################################
# Make a list of candidate genes to begin the prioritisation process

# Pull out all genes that overlap the locus boundary for all loci

loci = map(loci, get_boundary)
gene_map = map(loci, ~ get_overlap_genes(.x, gencode))

# Clean up eQTL data
eqtl_list = map(c('full', 'male', 'female'), ~ eqtl %>% filter(cohort == .x))
eqtl_list = map2(eqtl_list, loci, ~ clean_eqtl(.x, .y, gencode))

kidney_eqtl_list = clean_eqtl(kid_sub, loci[[1]], gencode)

# Add kidney eQTL to the full cohort eQTL result
eqtl_list[[1]] = rbind(eqtl_list[[1]], kidney_eqtl_list)

gene_map = map2(gene_map, eqtl_list, ~ add_eqtl_gene(.x, .y))

# Number of unique genes going into gene prioritisation:
full_list = map_dfr(gene_map, ~ .x)

full_list %>% filter(within_locus) %>% pull(gene) %>% unique %>% length
full_list %>% filter(eqtl != 'FALSE') %>% pull(gene) %>% unique %>% length
length(unique(full_list$gene))

################################################################################
# Pull out blood eQTLs
blood_eqtl = eqtl %>% filter(tissue == 'Whole_Blood') %>% select(SNP, gene_name, PP.H4.abf, cohort) %>% mutate(gtex_eqtl.Whole_Blood = T)
blood_eqtl_list = map(c('full', 'male', 'female'), ~ blood_eqtl %>% filter(cohort == .x) %>% select(!cohort))

# Add this to the gene list
gene_map = map2(gene_map, blood_eqtl_list, ~ left_join(.x, .y, by = c('SNP', 'gene' = 'gene_name')))

gene_map = map(gene_map, ~ .x %>% mutate( gtex_eqtl.Whole_Blood = ifelse(is.na(gtex_eqtl.Whole_Blood), F, gtex_eqtl.Whole_Blood), PP.H4.abf = ifelse(is.na(PP.H4.abf), "<0.8", PP.H4.abf)))

################################################################################
# Pull out those involved in meQTL

meqtl_list = map(c('full', 'male', 'female'), ~ meqtl %>% filter(cohort == .x) %>% group_by(SNP) %>% summarise(num_meQTL = n(), max_meqtl_PP.H4.abf = max(PP.H4.abf)) %>% ungroup %>% mutate(meqtl_locus = T))

# Add this to the gene list
gene_map = map2(gene_map, meqtl_list, ~ left_join(.x, .y, by = c('SNP')))

gene_map = map(gene_map, ~ .x %>% mutate( meqtl_locus = ifelse(is.na(meqtl_locus), F, meqtl_locus), num_meQTL = ifelse(is.na(num_meQTL), 0, num_meQTL), max_meqtl_PP.H4.abf = ifelse(is.na(max_meqtl_PP.H4.abf), "<0.8", max_meqtl_PP.H4.abf)))

################################################################################
# Add Ruth's pheWAS

phewas = read.table('/Volumes/archive/merrimanlab/major_gwas_paper_archive/results_for_paper/pheWAS/phewas_wbc.txt', sep = '\t', header = T, stringsAsFactors = F)
phewas = phewas %>% mutate(wbc_locus = T) %>% select(rsid, wbc_locus) %>% distinct

gene_map = map(gene_map, ~ left_join(.x, phewas, by = c('SNP' = 'rsid')))

gene_map = map(gene_map, ~ .x %>% mutate(wbc_locus = ifelse(is.na(wbc_locus), F, wbc_locus)))

################################################################################
# Differentially expressed genes in gout
degs = read.table('/Volumes/archive/merrimanlab/major_gwas_paper_archive/results_for_paper/misc/TopHits4Riku.csv', sep = ',', header = T, stringsAsFactors = F)

# Add gene positions
degs = gencode %>% select(Chrom:End, Gene, ensemblGeneID_short) %>% left_join(degs, ., by = c('EnsID' = 'ensemblGeneID_short'))

# CD24 is in chr6:107417706-107423502, according to ensembl
degs$Chrom[which(degs$GENEids == 'CD24')] = 'chr6'
degs$Start[which(degs$GENEids == 'CD24')] = 107417706
degs$End[which(degs$GENEids == 'CD24')] = 107423502
degs$Gene[which(degs$GENEids == 'CD24')] = 'CD24'

degs = degs %>% mutate(is_deg = T) %>% select(Gene, log2FoldChange, is_deg)

# Add this to the gene list
gene_map = map(gene_map, ~ left_join(.x, degs, by = c('gene' = 'Gene')))
gene_map = map(gene_map, ~ .x %>% mutate(is_deg = ifelse(is.na(is_deg), F, is_deg)))

################################################################################
# Gene has an eQTL (for the same EUR GWAS variant) in monocytes in either
# ImmuNexUT or OneK1K cohort

immun_ut = read.table('/Volumes/archive/merrimanlab/major_gwas_paper_archive/results_for_paper/misc/immunexut_eqtl.txt', sep = '\t', header = T, stringsAsFactors = F) %>% mutate(ENSG_short = gsub('\\..*', '', Gene_id))
immun_ut = immun_ut %>% select(Variant_ID, Gene_name, Forward_slope, ENSG_short) %>% mutate(immun_ut = T) %>% distinct
colnames(immun_ut)[1:3] = c('SNP', 'gene', 'slope')

onek1k = read.table('/Volumes/archive/merrimanlab/major_gwas_paper_archive/results_for_paper/misc/onek1k_all_snp_lead.fdr0.05.txt', sep = '\t', header = T, stringsAsFactors = F) %>% filter(grepl('Mono', CELL_TYPE)) %>% select(RSID, GENE, GENE_ID, SPEARMANS_RHO) %>% mutate(onek1k = T)

mono_eqtl = full_join(immun_ut, onek1k, by = c('SNP' = 'RSID', 'gene' = 'GENE', 'slope' = 'SPEARMANS_RHO', 'ENSG_short' = 'GENE_ID'))

mono_eqtl$immun_ut[is.na(mono_eqtl$immun_ut)] = F
mono_eqtl$onek1k[is.na(mono_eqtl$onek1k)] = F

mono_eqtl$monocyte_eqtl = mono_eqtl$immun_ut | mono_eqtl$onek1k

# Need to figure out which cohrot the variant came from:
loci_cohort = map2_dfr(loci, c('full', 'male', 'female'), ~ .x %>% mutate(cohort = .y) %>% distinct(SNP, cohort))

mono_eqtl = left_join(mono_eqtl, loci_cohort)

mono_eqtl_list = map(c('full', 'male', 'female'), ~ mono_eqtl %>% filter(cohort == .x) %>% select(!cohort) %>% group_by(SNP, gene) %>% mutate(max.slope = ifelse(sum(slope) < 0, min(slope), max(slope))) %>% select(!slope) %>% ungroup)

# Add this to gene list
gene_map = map2(gene_map, mono_eqtl_list, ~ .y %>% select(SNP, ENSG_short, max.slope, monocyte_eqtl) %>% left_join(.x, ., by = c('SNP', 'ENSG_short')))
gene_map = map2(gene_map, mono_eqtl_list, ~ .y %>% select(SNP, gene, max.slope, monocyte_eqtl) %>% left_join(.x, ., by = c('SNP', 'gene')))

gene_map = map(gene_map, ~ .x %>% mutate(monocyte_eqtl = monocyte_eqtl.x | monocyte_eqtl.y, max.slope = ifelse(is.na(max.slope.x), max.slope.y, max.slope.x)) %>% select(locus:is_deg, max.slope, monocyte_eqtl))

gene_map = map(gene_map, ~ .x %>% mutate(monocyte_eqtl = ifelse(is.na(monocyte_eqtl), F, monocyte_eqtl)) %>% distinct)

################################################################################
# Genes that are expressed in GTEx Whole Blood

gene_exp = read.table('/Volumes/archive/merrimanlab/major_gwas_paper_archive/results_for_paper/misc/gtex_whole_blood.gene_exp.txt', sep = '\t', header = T, stringsAsFactors = F)
gene_exp = gene_exp %>% select(Name:Description, mean_tpm, tpm_0.5)

colnames(gene_exp) = c("ENSG", "gene", "mean_tpm", "tpm_0.5")

gene_exp$ENSG_short = gsub('\\..*', '', gene_exp$ENSG)

# Match up with ENSG:
gene_map = map(gene_map, ~ left_join(.x, gene_exp, by = c('ENSG_short'), suffix = c('', '.y')) %>% select(!gene.y))

gene_map = map(gene_map, ~ .x %>% mutate(tpm_0.5 = ifelse(is.na(tpm_0.5), F, tpm_0.5)))

################################################################################
# Add DEG from MSU/LPS stimulated cells from Cobo et al.

cobo = read.table('/Volumes/archive/merrimanlab/major_gwas_paper_archive/results_for_paper/misc/cobo_deglist.combined.txt', sep = '\t', header = T, stringsAsFactors = F)

cobo = cobo %>% mutate(deg_msu_lps = T) %>% select(gene, contains('log2'), deg_msu_lps)

gene_map = map(gene_map, ~ left_join(.x, cobo))

gene_map = map(gene_map, ~ .x %>% mutate(deg_msu_lps = ifelse(is.na(deg_msu_lps), F, deg_msu_lps)))

################################################################################
# Add "Function-agnostic" score
################################################################################
# Add closest TSS score

tss = read.table('/Volumes/archive/merrimanlab/major_gwas_paper_archive/results_for_paper/misc/fantom5_tss.txt', sep = '\t', header = T, stringsAsFactors = F)

gene_list = map_dfr(gene_map, ~ .x %>% select(locus:snp_pos, ENSG_short)) %>% distinct

tss_map = left_join(gene_list, tss, by = c('snp_chr' = 'chr', 'ENSG_short' = 'ensemblGeneID')) %>% mutate(within_tss = (snp_pos >= tss_start & snp_pos <= tss_end), tss_dist = abs(snp_pos - tss_start))

closest_tss = tss_map %>% group_by(SNP) %>% arrange(tss_dist) %>% slice(1) %>% ungroup %>% mutate(closest_tss = !is.na(tss_start)) %>% select(locus:SNP, gene, tss_start:tss_end, tss_dist:closest_tss)

# Need to add TSS information first, then add info on whether the TSS is the
# closest for that particular SNP/locus

gene_map = map(gene_map, ~ tss_map %>% select(!gene) %>% left_join(.x, .))
gene_map = map(gene_map, ~ closest_tss %>% select(!gene) %>% left_join(.x, .))

gene_map = map(gene_map, ~ .x %>% mutate( within_tss = ifelse(is.na(within_tss), F, within_tss), closest_tss = ifelse(is.na(closest_tss), F, closest_tss)))

# Flag genes that contain lead variant
contain_lead = map(gene_map, ~ pmap_lgl(list(.x$snp_pos, .x$gene_start, .x$gene_end), ~ between(..1, ..2, ..3)))

gene_map = map2(gene_map, contain_lead, ~ .x %>% mutate(contain_lead = .y))

# Gene contains lead variant OR TSS is closest to lead variant (none of the
# lead SNPs were within TSS)
gene_map = map(gene_map, ~ .x %>% mutate(has_lead_or_tss = closest_tss | contain_lead))

################################################################################
# Get genes with evidence of regulation/connection by ABC-enhancers:

eqtl_abc = read.table('/Volumes/archive/merrimanlab/major_gwas_paper_archive/results_for_paper/abc/abc_eqtl_overlap.14JUL2023.eqtl_gene_match.txt', sep = '\t', header = T, stringsAsFactors = F)
eqtl_abc = left_join(eqtl_abc, loci_cohort)

nc_abc = read.table('/Volumes/archive/merrimanlab/major_gwas_paper_archive/results_for_paper/misc//non_coding_abc.txt', sep = '\t', header = T, stringsAsFactors = F)
nc_abc = nc_abc %>% separate_rows(cohort, sep = '/') %>% distinct

eqtl_abc = map(c('full', 'male', 'female'), ~ eqtl_abc %>% filter(cohort == .x) %>% select(TargetGene) %>% mutate(is_abc_target = T))
nc_abc = map(c('full', 'male', 'female'), ~ nc_abc %>% filter(cohort == .x) %>% select(TargetGene) %>% mutate(is_abc_target = T))

abc_list = map2(eqtl_abc, nc_abc, ~ rbind(.x, .y) %>% distinct)

gene_map = map2(gene_map, abc_list, ~ left_join(.x, .y, by = c('gene' = 'TargetGene')))
gene_map = map(gene_map, ~ .x %>% mutate(is_abc_target = ifelse(is.na(is_abc_target), F, is_abc_target)))

################################################################################
# Flag genes with missense variants

missense_genes = c( "ABCG2", "ADH1B", "ALDH2", "CPS1", "GCKR", "HNF4A", "MC4R", "PNPLA3", "SLC17A1", "SLC17A3", "SLC2A9", "ABCA6", "ADO", "AP4E1", "AQP10", "BSCL2", "CRIP3", "CUBN", "CTAGE9", "DTL", "DDIT4L", "EPB41", "EVI5", "FAM35A", "FGF21", "FRK", "GLIS3", "GLS2", "INHBC", "HNF1A", "JMJD1C", "KIAA0100", "LRP2", "MLXIPL", "MFSD12", "NPHS2", "POM121", "SH2B1", "SH2B3", "SLC25A5", "SLC5A1", "SLC5A9", "SLC39A8", "SLCO1B1", "SOS2", "TSPAN6", "UPF3A" )

missense_dat = data.frame(gene = missense_genes, candidate_missense = T)

gene_map = map(gene_map, ~ left_join(.x, missense_dat))
gene_map = map(gene_map, ~ .x %>% mutate(candidate_missense = ifelse(is.na(candidate_missense), F, candidate_missense)))

################################################################################
# Generate a prioritisation score

res = map(gene_map, ~ .x %>% select(locus:gene_end, tss_start:tss_end, tss_dist, ENSG_short:eqtl, gtex_eqtl.Whole_Blood, meqtl_locus:wbc_locus, is_deg, monocyte_eqtl, tpm_0.5, deg_msu_lps, PP.H4.abf, max_meqtl_PP.H4.abf, log2FoldChange, max.slope, mean_tpm, log2FoldChange.LPS:log2FoldChange.LPS_MSU, within_tss, closest_tss, is_abc_target, candidate_missense))

cols = c( "locus", "Main_GWAS_loci", "ancestry", "SNP", "snp_chr", "snp_pos", "gene", "gene_chr", "gene_start", "gene_end", "tss_start", "tss_end", "tss_dist", "ENSG_short", "coding", "within_locus", "eqtl", "gtex_eqtl.Whole_Blood", "meqtl_locus", "wbc_locus", "is_deg", "monocyte_eqtl", "gtex_exp.Whole_Blood", "deg_msu_lps", "blood_eqtl_PP.H4.abf", "max_meqtl_PP.H4.abf", "log2FoldChange.DEG", "max_slope.monocyte_eQTL", "mean_tpm.Whole_Blood", "log2FoldChange.LPS", "log2FoldChange.MSU", "log2FoldChange.LPS_MSU", "within_tss", "closest_tss", "is_abc_target", "candidate_missense" )

colnames(res[[1]]) = cols
colnames(res[[2]]) = cols
colnames(res[[3]]) = cols

prioritisation_score = map(res, ~ .x %>% select(gtex_eqtl.Whole_Blood:deg_msu_lps) %>% apply(., 1, function(x) sum(x, na.rm = T)))
fn_agnostic_score = map(res, ~ .x %>% select(closest_tss:candidate_missense) %>% apply(., 1, function(x) sum(x, na.rm = T)))

res = pmap(list(res, prioritisation_score, fn_agnostic_score), ~ ..1 %>% mutate(prioritisation_score = ..2, fn_agnostic_score = ..3))

res = map(res, ~ .x %>% arrange(desc(prioritisation_score), desc(fn_agnostic_score), gene_chr, gene_start) %>% distinct)

map(res, ~ table(.x$prioritisation_score))
map(res, ~ table(.x$fn_agnostic_score))

res = map2_dfr(res, names(res), ~ .x %>% mutate(cohort = .y) %>% select(locus:ancestry, cohort, SNP:cohort))

write.table(res, 'gene_prioritisation.txt', sep = '\t', col.names = T, row.names = F, quote = F)

################################################################################

# Normalise the scores based on how many categories are relevant for each
# SNP-gene pair

res_norm = res

# Extra notes on which categories should be counted for which ancestry:
# - if tss_start is not available, don't count tss
# - if EAS/AFR/LAT, don't count gtex_eqtl.Whole_Blood
# - don't count meQTL for non-EUR
#
# So, the final number of categories used for each ancestry will be:
# - EUR  = 7
# - AFR  = 5
# - EAS  = 5
# - LAT  = 5
# - TAMA = 6
#
# And the final number of function-agnostic categories used will be 3, but only
# for those genes that had TSS information. If no TSS info, then it should be 2.

col_eur = c("gtex_eqtl.Whole_Blood", "meqtl_locus", "wbc_locus", "is_deg", "monocyte_eqtl", "gtex_exp.Whole_Blood", "deg_msu_lps")

col_non_eur = c("wbc_locus", "is_deg", "monocyte_eqtl", "gtex_exp.Whole_Blood", "deg_msu_lps")

col_tama = c("gtex_eqtl.Whole_Blood", "wbc_locus", "is_deg", "monocyte_eqtl", "gtex_exp.Whole_Blood", "deg_msu_lps")

col_no_tss = c("is_abc_target", "candidate_missense")

# Re-score the prioritisation and function-agnnostic scores based on the new
# definition
eur_ind = which(res_norm$ancestry == 'EUR')
non_eur_ind = which(res_norm$ancestry %in% c('AFR', 'EAS', 'LAT'))
tama_ind = which(res_norm$ancestry == 'TAMA')

res_norm$prioritisation_score[eur_ind] = res_norm[eur_ind, ] %>% select(all_of(col_eur)) %>% apply(., 1, function(x) sum(x, na.rm = T))
res_norm$prioritisation_score[non_eur_ind] = res_norm[non_eur_ind, ] %>% select(all_of(col_non_eur)) %>% apply(., 1, function(x) sum(x, na.rm = T))
res_norm$prioritisation_score[tama_ind] = res_norm[tama_ind, ] %>% select(all_of(col_tama)) %>% apply(., 1, function(x) sum(x, na.rm = T))

# Now re-score function-agnostic score:
no_tss_ind = which(is.na(res_norm$tss_start))

res_norm$fn_agnostic_score[no_tss_ind] = res_norm[no_tss_ind, ] %>% select(all_of(col_no_tss)) %>% apply(., 1, function(x) sum(x, na.rm = T))

# Make a column with possible max-score:
res_norm$max_prioritisation_score = 7
res_norm$max_prioritisation_score[non_eur_ind] = 5
res_norm$max_prioritisation_score[tama_ind] = 6

res_norm$max_fn_agnostic_score = ifelse(is.na(res_norm$tss_start), 2, 3)

# Normalise the scores_norm based on max possible score:
res_norm$prioritisation_score_norm = (res_norm$prioritisation_score / res_norm$max_prioritisation_score) * 7
res_norm$fn_agnostic_score_norm = (res_norm$fn_agnostic_score / res_norm$max_fn_agnostic_score) * 3

res_norm = res_norm %>% arrange(desc(prioritisation_score_norm), desc(fn_agnostic_score_norm), snp_chr, snp_pos)

# There are some genes that came from eQTL data that doesn't have
# "Main_GWAS_loci" value
res_norm$Main_GWAS_loci[which(is.na(res_norm$Main_GWAS_loci))] = res_norm$locus[which(is.na(res_norm$Main_GWAS_loci))]

write.table(res_norm, 'gene_prioritisation.norm_score.txt', sep = '\t', col.names = T, row.names = F, quote = F)

# Make a unique gene list for the full table and per-ancestry
res_uniq = res_norm %>% group_by(gene) %>% arrange(desc(prioritisation_score_norm), desc(fn_agnostic_score_norm)) %>% slice(1) %>% ungroup %>% arrange(desc(prioritisation_score_norm), desc(fn_agnostic_score_norm), snp_chr, snp_pos)

write.table(res_uniq, 'gene_prioritisation.norm_score.unique_overall.txt', sep = '\t', col.names = T, row.names = F, quote = F)

