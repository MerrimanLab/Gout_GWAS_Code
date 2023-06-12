#! /bin/bash

# Script to generate baseline v1.1 LD score regression annotation to be used
# for 1KGP reference LD score regression

ANNOT_DIR=data/ldsc/ld_annot/baseline_v1.1
ANCESTRY=$1
OUT_DIR=data/ldsc/1kgp_ref_ldsc/${ANCESTRY}

for i in {1..22} ; do
	parallel -j 20 "make_annot.py --bed-file {} --bimfile data/ldsc/1kgp_plink/${ANCESTRY}/${ANCESTRY}_chr${i}.mac5.bim --annot-file ${OUT_DIR}/chr${i}.{/.}.tmp && sed 's:ANNOT:{/.}:g' ${OUT_DIR}/chr${i}.{/.}.tmp > ${OUT_DIR}/chr${i}.{/.}" ::: $(ls ${ANNOT_DIR}/*)
	rm ${OUT_DIR}/chr${i}.*.tmp
	awk 'BEGIN {print "CHR", "BP", "SNP", "CM", "base"}; {print $1, $4, $2, $3, 1}' data/ldsc/1kgp_plink/${ANCESTRY}/${ANCESTRY}_chr${i}.mac5.bim > ${OUT_DIR}/chr${i}_tmpbim
	paste ${OUT_DIR}/chr${i}_tmpbim ${OUT_DIR}/chr${i}.* | tr -s ' ' '\t' > ${OUT_DIR}/${ANCESTRY}_chr${i}.annot
	gzip -f ${OUT_DIR}/${ANCESTRY}_chr${i}.annot
	rm ${OUT_DIR}/chr${i}*
done

