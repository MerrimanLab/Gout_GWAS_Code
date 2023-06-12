#! /bin/bash

# Quick script to generate bim files to tack-on to the LD score annotations for
# CTG and CTS

KGP_DIR=data/ldsc/1kgp_plink/
CTG_DIR=data/ldsc/1kgp_ref_ctg/

for ancestry in {AFR,EAS,EUR,LAT} ; do
	for ct in {cts,ctg} ; do
		echo CHR BP SNP CM > ${CTG_DIR}/${ancestry}/${ct}/${ancestry}_head.txt
		for chr in {1..22} ; do
			cut -f1-4 ${KGP_DIR}/${ancestry}/${ancestry}_chr${chr}.mac5.bim | awk '{print $1, $4, $2, $3}' | cat ${CTG_DIR}/${ancestry}/${ct}/${ancestry}_head.txt - | tr ' ' '\t' > ${CTG_DIR}/${ancestry}/${ct}/${ancestry}_chr${chr}.bim
		done
	done
done

