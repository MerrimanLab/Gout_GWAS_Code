#! /bin/bash

# Script to calculate the 1000 Genomes MAF for LD score regression
# MAF-adjusting and binning

ANCESTRY=$1
BASE_DATA=data/ldsc/1kgp_plink/${ANCESTRY}/${ANCESTRY}_chr

parallel "plink1.9b4.9 --bfile ${BASE_DATA}{}.mac5 --freq --out ${BASE_DATA}{}.mac5" ::: {1..22}

