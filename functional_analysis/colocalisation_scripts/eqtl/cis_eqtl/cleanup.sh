# Remove files for lead snps that didn't have a significant eQTL
parallel -j 1 'rm {}_*.txt' ::: $(cat MISSING.non_significant.txt)

# Remove intermediate files
rm *subset*.txt *gene_lookup*.txt gwas_gtex.txt snplist.txt gtex_tissue_genes.txt gwas_sorted.txt *allpairs.txt*
