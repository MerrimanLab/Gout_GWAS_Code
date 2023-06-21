library(tidyverse)

args = commandArgs(trailingOnly=TRUE)

gwas_w_gtexids_file <- args[1]   # tab delimited gwas file with following columns: GTEX_ID CHR BP SNP N minor major MAF effect SE P
snp_list_file <- args[2]         # snp list. One snp per line. no header.
window_bp <- as.numeric(args[3]) # number - size of total window centered on lead snp

snp_list <- read_tsv(snp_list_file, col_names = c("snp"))

gwas <- vroom::vroom(gwas_w_gtexids_file, delim = "\t", col_types = cols( GTEX_ID = col_character(), CHR = col_double(),  BP = col_double(),  SNP = col_character(), N = col_double(), minor = col_character(),  major = col_character(),  MAF = col_double(), effect = col_double(),  SE = col_double(),  P = col_double()))

extract_region_snps_from_gwas <- function(snp, gwas, region_size_bp){
  half_window <- region_size_bp / 2
  snp_loc <- gwas %>% filter(SNP == snp) %>% select(CHR, BP)
  if(NROW(snp_loc) != 1){
    stop(paste0(snp," either wasn't found or matched mulitple rows in `gwas`"), call. = FALSE)
  }

  gwas %>% filter(snp_loc$CHR[[1]] == CHR & between(BP, snp_loc$BP[[1]] - half_window, snp_loc$BP[[1]] + half_window)) %>% arrange(CHR,BP) %>% select(GTEX_ID, everything()) %>%
    write_tsv(paste0(snp, "_gwas_snps.txt"))
}

safe_extract_region_snps_from_gwas <- purrr::safely(extract_region_snps_from_gwas)

purrr::walk(snp_list$snp, safe_extract_region_snps_from_gwas, gwas, window_bp)
