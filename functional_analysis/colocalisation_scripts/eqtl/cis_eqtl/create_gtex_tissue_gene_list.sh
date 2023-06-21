file=$1
rsid=$(basename -s _gwas_snps.txt ${file})

cut -f1 ${file} | grep -h -Fwf - gtex_tissue_genes.txt > ${rsid}_gtex_tissue_genes.txt

