#! /bin/bash

# Pull out variants from GWAS and meQTL data and merge it for colocalisation
# analysis
#
# Usage:
# bash src/prep_coloc_file.sh <path_to>/4_89052323_sig_mQTL.txt
#
# Note:
# - Merging of the two files (GWAS and mQTL data) will be done in the coloc script

CHR=$(echo ${1##*/} | sed 's/_.*//g')
SEX=$2
MQTL_PATH=/Volumes/archive/merrimanlab/reference_files/GoDMC_meQTL/cleaned/cis
GWAS=/Volumes/archive/merrimanlab/major_gwas_paper_archive/results_meta_analysis/EUR/${SEX}/EUR_meta_${SEX}1_clean_rsid.nfiltered.biallelic.txt

# Pull out all the variants that affect the mQTL sites listed in the input file
parallel -j15 "cat <(head -1 ${MQTL_PATH}/GoDMC_cis_meQTL.chr${CHR}.tsv) <(grep -Fw {} ${MQTL_PATH}/GoDMC_cis_meQTL.chr${CHR}.tsv) | grep -v 'INDEL' | cut -f1-4,6-13 | tee ${1%_sig_mQTL.txt}_{}.mqtl.txt | cut -f4 > ${1%_sig_mQTL.txt}_{}.lookup " ::: $(cat $1)

# Now pull out relevant variants from the GWAS data
parallel -j15 "grep -Fwf ${1%_sig_mQTL.txt}_{}.lookup ${GWAS} | cut -f1-6,10-12,18-19 > ${1%_sig_mQTL.txt}_{}.gwas.txt " ::: $(cat $1)

rm ${1%_sig_mQTL.txt}_cg*.lookup
