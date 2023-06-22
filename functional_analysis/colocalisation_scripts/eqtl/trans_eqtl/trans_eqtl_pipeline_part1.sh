################################################################################
#
# Pipeline to generate a list of trans-eQTLs from a SNP list and a gwas.
#
# bash trans_eqtl_pipeline_part1.sh <gwas> <qtl_list> <window size (bp)>
# bash trans_eqtl_pipeline_part1.sh my_gwas.txt my_qtl_list.txt 1000000
#
################################################################################

GWAS=$1      # gwas results file with columns according to the following header: CHR, BP, SNP, A1, A2, MAF, BETA, SE, P, N
QTL_LIST=$2  # tab delimited, header: RSID\tGENE\tTISSUE. Gene can be either ensembl gene id or HGNC symbol. As defined in gtex e.g. "Brain_Hippocampus"
WINDOW_BP=$3 # number for the window size (bp) to grab snps from surrounding lead snp e.g. 1000000 will look +500kb and -500kb centered on lead snp

# Test to ensure gwas header is correct
diff <(head -n 1 $GWAS | tr '\t' ' ') <(echo CHR BP SNP A1 A2 MAF BETA SE P N) > /dev/null
if [ $? -ne 0 ]
then
    echo "Header of $GWAS needs to match 'CHR BP SNP A1 A2 MAF BETA SE P N' (tab or space delimited)"
    echo "Exiting"
    exit
fi

################################################################################
# STEP 1: setup the qtl list and gwas files
################################################################################
echo "STEP 1 $(date) (about 5 min)"
diff <(head -n 1 $QTL_LIST) <(echo -e "RSID\tGENE\tTISSUE") > /dev/null
if [ $? -ne 0 ]
then
    echo "Header of $QTL_LIST needs to match 'RSID\tGENE\tTISSUE' (tab delimited)"
    echo "Exiting"
    exit
fi

# Take only rsid, tissue, and gene columns. Make sure only has unix line
# breaks, and no duplicates.
cut -f1-3 $QTL_LIST | tr -d '\r'| sort -u > qtl_list.txt

# The following command was how gtex_lookup_sorted_rsid.txt was originally created
# zcat GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.lookup_table.txt.gz | cut -f1,7 | sort -k 2b,2 -S 4G > gtex_lookup_sorted_rsid.txt

# Sort the gwas ready for joining
sort -k 3b,3 -S 8G ${GWAS} > gwas_sorted.txt

# Add the gtex id onto all the gwas snps where possible
join -1 2 -2 3 gtex_lookup_sorted_rsid.txt gwas_sorted.txt 2> gwas_gtex.log | cat <(echo "SNP GTEX_ID CHR BP minor major MAF effect SE P N") -  | tr ' ' '\t' > gwas_gtex.txt
# gwas_gtex.txt has header of (tab delimited):
# SNP, GTEX_ID, CHR, BP, N, A1, A2, MAF, BETA, SE, P

################################################################################
# STEP 2: Find all the SNPs from the GWAS we want
################################################################################

echo "STEP 2 $(date) (about 30sec / snp)"

# From the list of snps, create a file per lead SNP that has all the other snps from the gwas within 1 megabase
# (takes a while - depends on number of lead and gwas snps, ~30sec/lead snp)
Rscript gwas_extract_snps_in_region.R gwas_gtex.txt qtl_list.txt $WINDOW_BP 2> surrounding_snps.log

################################################################################
# STEP 3 Work out what lead snps aren't going to be used and why
################################################################################

# Create list of lead snps that were in GTEx
echo "STEP 3 $(date)"
cut -f1 qtl_list.txt | sed '1d' | sort -u > snplist.tqtl.txt
grep -Fwf snplist.tqtl.txt gtex_lookup_sorted_rsid.txt > tmp_snplist_matched.tqtl.txt
cut -f1 tmp_snplist_matched.tqtl.txt > tmp_gtexid.tqtl.txt
cut -f2 tmp_snplist_matched.tqtl.txt > tmp_rsid.tqtl.txt

# Use this list to find the snps that weren't in GTEx
grep -v -Fwf tmp_rsid.tqtl.txt snplist.tqtl.txt > MISSING.no_gtex_match.tql.txt

# Create a list of snps to be removed because the eQTL is non-sigificant
# Create list of SNPs that weren't in the GWAS
tr ' ' '\t' < ${GWAS} | cut -f 3 > tmp_gwas_rsids.tqtl.txt
grep -v -Fwf tmp_gwas_rsids.tqtl.txt snplist.tqtl.txt | grep -v "^RSID"> MISSING.gwas.tqtl.txt

# Remove files for lead snps that didn't have a significant eQTL
rm tmp_*

if [  $(ls rs*.gwas.txt | wc -l) -eq 0 ]
then
    echo All lead SNPs were either missing from GTex or non-sigificant
    echo check MISSING.non_significant.tqtl.txt and MISSING.no_gtex_match.tqtl.txt
    echo EXITING
    exit
fi

echo "Starting list had $(wc -l snplist.tqtl.txt) lead SNPs (after removing duplicates)"
echo "$(grep -Fwf snplist.tqtl.txt ${GWAS} | wc -l) lead SNPs matched with the GWAS file"
echo "$(wc -l MISSING.gwas.tqtl.txt) lead SNPs didn't match with the GWAS (MISSING.gwas.tqtl.txt)"
echo "$(wc -l MISSING.no_gtex_match.tqtl.txt) lead SNPS had no match in GTEx so were left out (MISSING.no_gtex_match.tqtl.txt)."

echo "End of part 1. Collect the rs*_trans_eqtl_snps.list.txt files and send them off. Once you have the results back start part 2."

