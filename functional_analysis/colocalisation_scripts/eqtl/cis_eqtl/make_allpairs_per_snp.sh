file=$1
rsid=$(basename -s _gwas_snps.txt ${file})

grep -h -Fwf <(cut -f1 ${file}) *.allpairs_subset.txt > ${rsid}.allpairs.txt
