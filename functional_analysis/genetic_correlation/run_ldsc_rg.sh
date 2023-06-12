#! /bin/bash

# Script to run LDSC genetic correlation with Neale UKBB phenotypes

TRAIT1=$1
TRAIT2=$2
TRAIT_BASE=${TRAIT2##*/}
OUT_DIR=$3

# Use "eur_w_ld_chr" for both ref and weights
ldsc.py --rg ${TRAIT1},${TRAIT2} --ref-ld-chr data/ldsc/ref_files/eur_w_ld_chr/ --w-ld-chr data/ldsc/ref_files/eur_w_ld_chr/ --chisq-max 80 --out ${OUT_DIR}/${TRAIT_BASE%%.*}
