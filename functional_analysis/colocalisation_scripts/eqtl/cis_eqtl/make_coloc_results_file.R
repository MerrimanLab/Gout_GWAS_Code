library(tidyverse)
library(snpStats)
library(coloc)

args = commandArgs(trailingOnly=TRUE)

file <- args[1]
lead_snp <- str_remove(file, "_gwas_snps.txt")
gtex_file <- paste0(lead_snp,".allpairs.txt")
sig_snp_gene_tissues <- read_delim("gtex_tissue_genes.txt", col_names = c("tissue", "GTEX_ID","gene_id"), delim =" ")
gwas_snps <- read_tsv(file)

# Establish which tissues the lead snp was significant in
lead_snp_sig_gene_tissues <- gwas_snps %>% select(SNP, GTEX_ID) %>% filter(SNP == lead_snp) %>% distinct() %>% inner_join(sig_snp_gene_tissues) %>% select(-SNP)

tissues <- read_tsv(gtex_file, col_names = c("tissue","gene_id", "variant_id", "tss_distance",  "ma_samples", "ma_count", "maf", "pval_nominal", "slope", "slope_se"),
                    col_types = cols(
                      tissue = col_character(),
                      gene_id = col_character(),
                      variant_id = col_character(),
                      tss_distance = col_double(),
                      ma_samples = col_double(),
                      ma_count = col_double(),
                      maf = col_double(),
                      pval_nominal = col_double(),
                      slope = col_double(),
                      slope_se = col_double()
                    )
    ) %>% inner_join(lead_snp_sig_gene_tissues, by = c("tissue", "gene_id"))

# Make a gwas file per gene per tissue
gwas_gtex <- left_join(gwas_snps, tissues, by = c("GTEX_ID"= "variant_id")) %>%
  select(GTEX_ID,
         tissue,
         gene_id,
         RSID=SNP,
         GWAS_Effect=effect,
         GWAS_SE = SE,
         GWAS_P_value = P,
         GWAS_n_total_sum = N,
         eQTL_pval_nominal = pval_nominal,    eQTL_ma = ma_count,
         eQTL_slope = slope,
         eQTL_slope_se = slope_se,
         eQTL_ma_samples = ma_samples,
         GWAS_maf = MAF,
         eQTL_maf = maf)  %>%
  mutate(eQTL_varbeta = (eQTL_slope_se) ^ 2,
         GWAS_varbeta = (GWAS_SE) ^ 2) %>%
  filter(eQTL_slope_se != "NaN")

write_tsv(gwas_gtex, paste0("Coloc_results/",lead_snp,".coloc_gwas_gtex_data.txt"))

# Establish which tissues the lead snp was significant in
combos <- lead_snp_sig_gene_tissues %>% select(tissue, gene_id) %>% distinct()

run_coloc <- function(tissue_name, gene_name, type = 'quant', prop_case = 0){
	coloc_data <- gwas_gtex %>% filter(tissue == tissue_name, gene_id == gene_name)

	# gwas data - change the data type, depending on quantitative or case
	# control data
	if (type == 'cc') {
		if (prop_case == 0) {
			stop('Invalid case proportion number')
		}
        message("Coloc using cc mode")
		gwas = list(snp = coloc_data$RSID,
					beta = coloc_data$GWAS_Effect,
					MAF = coloc_data$GWAS_maf,
					N = coloc_data$GWAS_n_total_sum,
					varbeta = coloc_data$GWAS_varbeta,
					pvalues = coloc_data$GWAS_P_value,
					type = type,
					s = prop_case)
	} else {
        message("Coloc using quant mode")
		gwas = list(snp = coloc_data$RSID,
					beta = coloc_data$GWAS_Effect,
					MAF = coloc_data$GWAS_maf,
					N = coloc_data$GWAS_n_total_sum,
					varbeta = coloc_data$GWAS_varbeta,
					pvalues = coloc_data$GWAS_P_value,
					type = type)
	}

	# gtex data
	gtex = list(snp = coloc_data$RSID,
				beta = coloc_data$eQTL_slope,
				MAF = coloc_data$eQTL_maf,
				N = coloc_data$eQTL_ma_samples,
				varbeta = coloc_data$eQTL_varbeta,
				pvalues = coloc_data$eQTL_pval_nominal,
				type = "quant")

	results <- coloc.abf(gwas, gtex)
	tibble(names = names(results$summary),values = results$summary) %>% mutate(gene_id = gene_name, tissue = tissue_name) %>% pivot_wider(tissue:gene_id, names_from = names, values_from = values)
}

study_type = args[2]
proportion = as.numeric(args[3])

#tictoc::tic()
results <- map2_dfr(combos$tissue, combos$gene_id, ~ run_coloc(.x, .y, type = study_type, prop_case = proportion))
#tictoc::toc()

# Save results
results %>% mutate(snp = lead_snp) %>% select(snp, everything())%>% write_tsv(paste0("Coloc_results/",lead_snp,".coloc_summary.txt"))
