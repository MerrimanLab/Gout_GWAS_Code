#! /bin/bash

# Script to make the CTG annotations for partitioned LD score
# regression

PREFIX=$1
ANCESTRY=$(echo ${PREFIX} | cut -d'_' -f1)
CTG=$(echo ${PREFIX} | cut -d'.' -f2)
CHR=$(echo ${PREFIX} | cut -d'.' -f3)
BED_FILE=data/ldsc/ld_annot/ct_and_ctg_bedfiles/cell_type_group_specific/${CTG}.bed
KGP_FILE=data/ldsc/1kgp_plink/${ANCESTRY}/${ANCESTRY}_chr${CHR}.mac5.bim
CTG_DIR=data/ldsc/1kgp_ref_ctg/${ANCESTRY}/ctg
BIM_FILE=${CTG_DIR}/${ANCESTRY}_chr${CHR}.bim
ANNOT_PRE=${CTG_DIR}/${ANCESTRY}_ctg.${CTG}.${CHR}

# Make cell type group annotation:
make_annot.py --bed-file ${BED_FILE} --bimfile ${KGP_FILE} --annot-file ${ANNOT_PRE}.tmp
paste ${BIM_FILE} ${ANNOT_PRE}.tmp > ${ANNOT_PRE}.annot
gzip -f ${ANNOT_PRE}.annot && rm ${ANNOT_PRE}.tmp

