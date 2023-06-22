library(tidyverse)
library(glue)
library(snpStats)
library(coloc)

args = commandArgs(trailingOnly=TRUE)

# Project directory
path <- paste0(args[1],"/")

message(paste0("Working in ",path))

# Takes a full file name and a project dir
#
# Needs the results of the tqtl to bein part2_input/ and have the same prefix
# names matching the gwas_results/ files
do_coloc <- function(full_filename, project_dir){
  filename <- basename(full_filename)

  prefix <- str_remove(filename, "-trans_eqtl_snps.list.txt")

  tqtl_file <- paste0(project_dir,"part2_input/",filename)
  gwas_file <- paste0(project_dir,"gwas_results/",prefix, "-trans_eqtl_snps.gwas.txt")

  # Check if both files exist, if not return NAs
  if(!(file.exists(tqtl_file) & file.exists(gwas_file))){
    return(tibble(names = c("nsnps", glue::glue("PP.H{h}.abf", h = 0:4)), values = NA, snp_gene_tissue = prefix) %>% pivot_wider(names_from = "names", values_from = values) )
  } else {

  transqtl_list <- read_tsv(tqtl_file)
  gwas_snps <- read_tsv(file = gwas_file)

  marker_list <- intersect(gwas_snps$GTEX_ID, transqtl_list$variant_id)
  gwas_snps <- gwas_snps %>% filter(GTEX_ID %in% marker_list)
  transqtl_list <- transqtl_list %>% filter(variant_id %in% marker_list)

  results <- coloc.abf(
    list(snp = gwas_snps$GTEX_ID,
         beta = gwas_snps$effect,
         MAF = gwas_snps$MAF,
         N = gwas_snps$N,
         varbeta = (gwas_snps$SE) ^ 2,
         pvalues = gwas_snps$P,
         type = "quant"),
    # gtex data
    list(snp = transqtl_list$variant_id,
         beta = transqtl_list$beta,
         MAF = transqtl_list$maf,
         N = transqtl_list$n,
         varbeta = (transqtl_list$beta_se) ^ 2,
         pvalues = transqtl_list$pval,
         type = "quant")
  )

  tibble(names = names(results$summary),values = results$summary, snp_gene_tissue = prefix) %>% pivot_wider(names_from = "names", values_from = values)  %>% write_tsv(paste0(project_dir,"Coloc_results/",prefix,".tqtl.coloc_summary.txt"))
  }
}

files <- list.files(path = paste0(path, "part2_input/"), pattern = "*-trans_eqtl_snps.list.txt", full.names = TRUE)
if(!dir.exists(paste0(path,"Coloc_results/"))){
  dir.create(paste0(path,"Coloc_results/"))
}

all_results <- map_dfr(files, ~do_coloc(.x, project_dir = path))
all_results %>% write_tsv(paste0(path,"Coloc_results/all_snps.tqtl.coloc_summary.txt"))

