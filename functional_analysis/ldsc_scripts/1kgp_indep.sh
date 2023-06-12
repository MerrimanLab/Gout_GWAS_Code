#! /bin/bash

# Script to make an LD-independent 1KGP data set that will be used for PC
# calculation.

ANCESTRY=$1
BASE_DATA=data/1kgp_ref/${ANCESTRY}_wgs
OUT_DIR=data/ldsc/1kgp_plink/${ANCESTRY}

# Make an LD-independent set of variants for PC calculation:
plink1.9b4.9 --bfile ${BASE_DATA} --indep-pairwise 50 5 0.2 --out ${OUT_DIR}/${ANCESTRY}_wgs_indep
plink1.9b4.9 --bfile ${BASE_DATA} --extract ${OUT_DIR}/${ANCESTRY}_wgs_indep.prune.in --make-bed --out ${OUT_DIR}/${ANCESTRY}_wgs_indep
plink1.9b4.9 --bfile ${OUT_DIR}/${ANCESTRY}_wgs_indep --recode --out ${OUT_DIR}/${ANCESTRY}_wgs_indep

