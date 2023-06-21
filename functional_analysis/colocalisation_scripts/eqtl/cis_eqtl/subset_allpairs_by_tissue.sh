tissue=$1

# Make a file per tissue of all the snps/genes we want
grep -w ${tissue} gtex_tissue_genes.txt | cut -d ' ' -f3 | sort -u > ${tissue}.gene_lookup.txt

# Pull out the relevant genes from the original data
zgrep -Fwf ${tissue}.gene_lookup.txt GTEx_eQTL_v8/Original_GTEx_files/${tissue}.allpairs.txt.gz | awk -v tissue=${tissue} '{print tissue"\t"$0}'> ${tissue}.allpairs_subset.txt
