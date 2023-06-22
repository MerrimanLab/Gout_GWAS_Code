echo "Collecting Coloc results"

cat Coloc_results/rs*summary.txt | grep -v "^snp" | sed 's/-/\t/' | sed 's/-/\t/' | cat <(echo -e  "snp\tgene_id\ttissue\tnsnps\tPP.H0.abf\tPP.H1.abf\tPP.H2.abf\tPP.H3.abf\tPP.H4.abf" ) - > Coloc_results/coloc_all_snps.txt

# Add gene names on
Rscript -e "library(dplyr); library(readr) ; read_tsv('Coloc_results/coloc_all_snps.txt') %>% left_join(read_tsv('GTEx_Analysis_v8_eQTL/ensembl_gene_lookup.txt'), by = c('gene_id' = 'gene_name')) %>% arrange(desc(PP.H4.abf)) %>% rename(ENSG = gene_id.y) %>% write_tsv('Coloc_results/all_snps_with_genes.tqtl.coloc_summary.txt')"

