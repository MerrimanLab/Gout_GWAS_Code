#! /bin/bash

# Script to generate reference LD scores for CTG using cov-LDSC

PREFIX=$1
ANCESTRY=$(echo ${PREFIX} | cut -d'_' -f1)
CHR=$(echo ${PREFIX} | cut -d'.' -f3)
KGP_FILE=data/ldsc/1kgp_plink/${ANCESTRY}/${ANCESTRY}_chr${CHR}.mac5
PC=data/ldsc/1kgp_pc/${ANCESTRY}_eigen.pca.evec.clean
CTG_DIR=data/ldsc/1kgp_ref_ctg/${ANCESTRY}/ctg
ANNOT=data/ldsc/1kgp_ref_ctg/${ANCESTRY}/ctg/${PREFIX}.annot.gz

time cov-ldsc.py --bfile ${KGP_FILE} --annot ${ANNOT} --cov ${PC} --print-snps data/ldsc/ref_files/w_hm3_nomhc_snplist_plink.txt --chunk-size 2000 --l2 --ld-wind-cm 20 --out ${CTG_DIR}/${PREFIX}

