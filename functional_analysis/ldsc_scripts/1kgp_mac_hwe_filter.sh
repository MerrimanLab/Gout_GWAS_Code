#! /bin/bash

# Script to split out the 1000 Genomes WGS PLINK file into chromosomes, and
# make a filtered and subsetted data set that will be used for LD score
# regression.

ANCESTRY=$1
BASE_DATA=data/1kgp_ref/${ANCESTRY}_wgs
OUT_DIR=data/ldsc/1kgp_plink/${ANCESTRY}

# Split it out into chromosomes, filter out MAC < 5 and put a HWE P < 1e-3
# threshold:
parallel -j 10 "plink1.9b4.9 --bfile ${BASE_DATA} --chr {} --mac 5 --hwe 1e-3 --cm-map data/ldsc/ref_files/genetic_map_chr{}_combined_b37.txt {} --make-bed --out ${OUT_DIR}/${ANCESTRY}_chr{}.mac5" ::: {{1..16},{18..22}}

# NOTE: There is one variant in chr 17 that is duplicated and cause problem
# when generating annotations in AFR and TAMA 1KGP data:
if [[ "${ANCESTRY}" == "AFR" || "${ANCESTRY}" == "TAMA" ]]; then
	plink1.9b4.9 --bfile ${BASE_DATA} --chr 17 --mac 5 --hwe 1e-3 --cm-map data/ldsc/ref_files/genetic_map_chr17_combined_b37.txt 17 --exclude-snp 17_1144632 --make-bed --out ${OUT_DIR}/${ANCESTRY}_chr17.mac5
else
	plink1.9b4.9 --bfile ${BASE_DATA} --chr 17 --mac 5 --hwe 1e-3 --cm-map data/ldsc/ref_files/genetic_map_chr17_combined_b37.txt 17 --make-bed --out ${OUT_DIR}/${ANCESTRY}_chr17.mac5
fi

# Copy the .frq file with a different file name.
# (It's slightly inconvenient to use the above naming for LDSC. Rather than
# renaming everything leading up to the LDSC analysis (e.g. making reference
# LDSC), I decided to copy the .frq file with a different name so I don't have
# to re-run everything - 15/1/2021)

for i in {1..22} ; do
	cp ${OUT_DIR}/${ANCESTRY}_chr${i}.mac5.frq ${OUT_DIR}/${ANCESTRY}_mac5_chr${i}.frq
done

