################################################################################

# Pipeline to generate a list of cis-eQTLs from a SNP list and a gwas run as
# following:

# quantitative trait:
# bash cis_eqtl_pipeline.sh <gwas> <snplist> <window size> <number of cpus> <study_type>
# bash cis_eqtl_pipeline.sh my_gwas.txt my_snplist.txt 1000000 8 quant 1> cis_pipeline.log 2>> cis_pipeline.log
#
# case-control:
# bash cis_eqtl_pipeline.sh <gwas> <snplist> <window size> <number of cpus> <study_type> <case_proportion>
# bash cis_eqtl_pipeline.sh my_gwas.txt my_snplist.txt 1000000 8 cc 0.05 1> cis_pipeline.log 2>> cis_pipeline.log

################################################################################

GWAS=$1       # gwas results file with columns according to the following header: CHR, BP, SNP, A1, A2, MAF, BETA, SE, P, N
SNPLIST=$2    # tab delimited, no header, snps in column 1
WINDOW_BP=$3  # number for the window size (bp) to grab snps from surrounding lead snp e.g. 1000000
CPUS=$4       # number of cpus for the parallel tasks - don't make too high because these are data read/write intensive
              # N.B. last stage uses $CPU for number of jobs, each job can also be assigned a number of cpus (default was 4) so min cpus was 32 for pipeline
STUDY_TYPE=$5 # Either 'quant' or 'cc'; if 'cc', then you must provide an extra argument (case proportion value)

# Check what the study type is. If it's case-control, then also check if it has
# a valid case proportion number:
if [[ ${STUDY_TYPE} == 'quant' ]]; then
    CASE_PROP=0
elif [[ ${STUDY_TYPE} == 'cc' ]]; then
    if [[ -z $6 ]]; then
        echo "You provided a case-control data without supplying a value for case_proportion - exiting"
        exit 1
    fi
    CASE_PROP=$6
fi

echo COMMAND USED: $0 $1 $2 $3 $4 $5 $6

# For meta results you can process with:
# awk '{print $1,$2,$3,$4,$5,$6,$10,$11,$12,$18}' < meta_gwas.txt > my_gwas.txt
# to make the correct order
#
# Original header: CHR BP SNP minor major MAF FreqSE MinFreq MaxFreq effect SE P Direction HetISq HetChiSq HetDf HetPVal N cpid SNP_original OR log10P

# Test to ensure gwas header is correct
diff <(head -n 1 $GWAS | tr '\t' ' ') <(echo CHR BP SNP A1 A2 MAF BETA SE P N) > /dev/null
if [ $? -ne 0 ]
then
    echo "Header of $GWAS needs to match 'CHR BP SNP A1 A2 MAF BETA SE P N' (tab or space delimited)"
    echo "Exiting"
    exit
fi

################################################################################
# STEP 1: setup the snp list and gwas files
################################################################################

echo "STEP 1 $(date) (about 5 min)"
cut -f1 $SNPLIST | tr -d '\r'| sort -u > snplist.txt

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
# From the list of snps, create a file per lead SNP that has all the other snps
# from the gwas within 1 megabase (takes a while - depends on number of lead
# and gwas snps, ~30sec/lead snp)
Rscript gwas_results_surrounding_snps.R gwas_gtex.txt snplist.txt $WINDOW_BP 2> surrounding_snps.log

################################################################################
# STEP 3: create lists of the relevant tissues/genes for each gwas list and remove lead snp files that won't/can't be persued
################################################################################

echo "STEP 3 $(date)"
# Combine all snps that are to be used and create a subset file containing
# tissue/gtex_id/gene so that this can be searched against to find the relevant
# gtex files to then go and search for our data. This is mostly for efficiency.
cat rs*_gwas_snps.txt | cut -f1 | sort -u | grep -h -Fwf - GTEx_Analysis_v8_eQTL/SNP_search_files/*txt | sort -u > gtex_tissue_genes.txt

# STEP 3a
# Work out what lead snps aren't going to be used and why

# Create list of lead snps that were in GTEx
echo "STEP 3a $(date)"
grep -Fwf snplist.txt gtex_lookup_sorted_rsid.txt > tmp_snplist_matched.txt
cut -f1 tmp_snplist_matched.txt > tmp_gtexid.txt
cut -f2 tmp_snplist_matched.txt > tmp_rsid.txt

# Use this list to find the snps that weren't in GTEx
grep -v -Fwf tmp_rsid.txt snplist.txt > MISSING.no_gtex_match.txt

# Create a list of snps to be removed because the eQTL is non-sigificant
grep -Fwf tmp_gtexid.txt gtex_tissue_genes.txt | cut -d' ' -f2 | sort -u > tmp_sig_gtexid.txt
grep -v -Fwf tmp_sig_gtexid.txt tmp_snplist_matched.txt | cut -f2 > MISSING.non_significant.txt

# Create list of SNPs that weren't in the GWAS
tr ' ' '\t' < ${GWAS} | cut -f 3 > tmp_gwas_rsids.txt
grep -v -Fwf tmp_gwas_rsids.txt snplist.txt > MISSING.gwas.txt

# Remove files for lead snps that didn't have a significant eQTL
rm tmp_*

if [  $(ls rs*_gwas_snps.txt | wc -l) -eq 0 ]
then
    echo All lead SNPs were either missing from GTex or non-sigificant
    echo check MISSING.non_significant.txt and MISSING.no_gtex_match.txt
    echo EXITING
    exit
fi

echo "Starting list had $(wc -l snplist.txt) lead SNPs (after removing duplicates)"
echo "$(grep -Fwf snplist.txt ${GWAS} | wc -l) lead SNPs matched with the GWAS file"
echo "$(wc -l MISSING.gwas.txt) lead SNPs didn't match with the GWAS (MISSING.gwas.txt)"
echo "$(wc -l MISSING.non_significant.txt) lead SNPs had non-significant eQTLs, and were removed (MISSING.non_significant.txt)"
echo "$(wc -l MISSING.no_gtex_match.txt) lead SNPS had no match in GTEx so were left out (MISSING.no_gtex_match.txt)."

# STEP 3b
# Search subsetted file to create a list of tissues and genes for each lead snp

echo STEP 3b $(date)
parallel -j $CPUS 'bash create_gtex_tissue_gene_list.sh {}' ::: rs*_gwas_snps.txt

################################################################################
# STEP 4: grab the info from the GTEx_eQTL_v8/ allpairs for the snps
################################################################################

echo "STEP 4 $(date) (~ 20 min)"

# Make a file per tissue of all the snps/genes we want
#
# Takes the gene column from the combined gtex_tissue_genes.txt file
# (takes about a 20min)
parallel -j $CPUS 'bash subset_allpairs_by_tissue.sh {}' ::: $(cut -d ' ' -f1 gtex_tissue_genes.txt | sort -u)

# Make a file per lead snp + gene + tissue combo and store them in directories
# by per lead snp inside Coloc_input/
#
# Input files are ${rsid}_gwas_snps.txt

# About 15 min / 8 snps
echo -e "Making snp allpairs files (Can take a couple of hours)"
parallel -j $CPUS 'bash make_allpairs_per_snp.sh {}' ::: rs*_gwas_snps.txt

################################################################################
# STEP 5: Run Coloc
################################################################################

mkdir -p Coloc_results
echo "STEP 5: Running coloc $(date) (~ 5min )"

# Run Coloc using $CPU jobs (default 8) each job using 4 cores
#
# ~25 min per snp
export OPENBLAS_NUM_THREADS=1
export STUDY_TYPE
export CASE_PROP
parallel -j $CPUS "Rscript make_coloc_results_file.R {} \$STUDY_TYPE \$CASE_PROP" ::: rs*_gwas_snps.txt

echo "Collecting Coloc results"
cat Coloc_results/rs*summary.txt | grep -v "^snp" | cat <(echo -e  "snp\ttissue\tgene_id\tnsnps\tPP.H0.abf\tPP.H1.abf\tPP.H2.abf\tPP.H3.abf\tPP.H4.abf" ) - > Coloc_results/coloc_all_snps.txt

# Add gene names on
Rscript -e "library(dplyr); library(readr) ; read_tsv('Coloc_results/coloc_all_snps.txt') %>% left_join(read_tsv('GTEx_Analysis_v8_eQTL/ensembl_gene_lookup.txt'), by = c('gene_id' = 'gene_id')) %>% arrange(desc(PP.H4.abf)) %>% write_tsv('Coloc_results/coloc_all_snps_with_genes.txt')"

echo "Combined results file 'Coloc_results/coloc_all_snps_with_genes.txt' made"

echo "Cleaning up"
bash cleanup.sh

echo "DONE $(date)"

