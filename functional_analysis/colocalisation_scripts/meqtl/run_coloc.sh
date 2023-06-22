#! /bin/bash

# Runs colocalisation analysis of GWAS with meQTL, given a file with list of
# significantly altered meQTL sites (from get_sig_meqtl.sh)
#
# Usage:
# bash src/run_coloc.sh 4_89052323_sig_mQTL.txt

export OMP_NUM_THREADS=1
export GOTO_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

module load bgenix/bgenix-1.1.8
module load bcftools/bcftools-1.11
module load plink/plink1.9b6.10

for i in $(cat $1); do
Rscript src/meqtl_scripts/coloc_script.R ${1%_sig_mQTL.txt}_${i}.gwas.txt ${1%_sig_mQTL.txt}_${i}.mqtl.txt
done

