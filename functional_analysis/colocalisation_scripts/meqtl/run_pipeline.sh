#! /bin/bash

export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OMP_NUM_THREADS=1
export GOTO_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

module load bgenix/bgenix-1.1.8
module load bcftools/bcftools-1.11
module load plink/plink1.9b6.10

# Given a list of CHR_BP, run the EUR gout GWAS vs GoDMC meQTL colocalisation
# analysis
#
# Usage:
# bash src/run_pipeline.sh <chr_bp_list> <full/male/female> <output_dir>

mkdir -p ${3}/{all_sig_mqtls,lz_plots}/${2}

# Pull out all lead variant-meQTL information from GoDMC data set:
parallel -j 100 "bash src/meqtl_scripts/get_sig_meqtl.sh {} ${3}/all_sig_mqtls/${2}/" ::: $(cat $1)

# Remove files with no lines, since files with no results cause problems when
# reading it in R (in the next command):
find ${3}/all_sig_mqtls/${2} -size 0 -print | grep -v 'not_found' > ${3}/all_sig_mqtls/${2}/not_found_cpid.txt
find ${3}/all_sig_mqtls/${2} -size 0 -delete

# Make sure the meQTLs are in cis (<= 1Mb from lead) and passes Bonferroni corrected p-value:
Rscript src/meqtl_scripts/check_sig_mqtl.R ${3}/all_sig_mqtls/${2}

parallel -j 15 "bash src/meqtl_scripts/prep_coloc_file.sh {} ${2}" ::: $(ls ${3}/all_sig_mqtls/${2}/*[0-9]_sig_mQTL.txt)

parallel -j 50 "bash src/meqtl_scripts/run_coloc.sh {}" ::: $(ls ${3}/all_sig_mqtls/${2}/*[0-9]_sig_mQTL.txt)

bash src/meqtl_scripts/combine_coloc_summary.sh ${3}/all_sig_mqtls/${2}/

# rm ${3}/*cg*
