#!/bin/bash

# Script to calculate LD score regression weights from 1KGP reference panel
ANCESTRY=$1
OUT_DIR=data/ldsc/1kgp_ref_ldsc/${ANCESTRY}

parallel "ldsc.py --bfile data/ldsc/1kgp_plink/${ANCESTRY}/${ANCESTRY}_chr{}.mac5 --l2 --ld-wind-cm 1 --extract data/ldsc/ref_files/w_hm3_nomhc_snplist_plink.txt --out ${OUT_DIR}/${ANCESTRY}_weights_chr{}" ::: {1..22}
