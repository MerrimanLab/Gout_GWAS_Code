#! /bin/bash

# A script to merge the X chromosome information onto the rest of the summary
# stat, from the output of SNPTEST

AUTO=${1}_autosomes.out

# Clean autosome data:
Rscript src/snptest_allele_correction.R ${AUTO}
awk 'NR > 1 {print $1, $3, $4, $5, $6, $18, $23, $28, $53, $54, $55, $45, $56, $48}' ${AUTO%.*}_corrected.txt > ${AUTO%.*}_clean.out

CHRX=${1}_chrX.out

# Clean chr X data and re-label X and/or XY to 23:
if [[ -f "$CHRX" ]]; then
	Rscript src/snptest_allele_correction_X.R ${CHRX}
	awk 'NR > 1 {print $1, "23", $4, $5, $6, $14, $24, $34, $50, $51, $52, $44, $53, $41 }' ${CHRX%.*}_corrected.txt > ${CHRX%.*}_clean.out
else
	PAR=${1}_PAR.out
	NONPAR=${1}_nonPAR.out
	Rscript src/snptest_allele_correction.R ${PAR}
	awk 'NR > 1 {print $1, "23", $4, $5, $6, $18, $23, $28, $53, $54, $55, $45, $56, $48}' ${PAR%.*}_corrected.txt > ${PAR%.*}_clean.out
	Rscript src/snptest_allele_correction_X.R ${NONPAR}
	awk 'NR > 1 {print $1, "23", $4, $5, $6, $14, $24, $34, $50, $51, $52, $44, $53, $41 }' ${NONPAR%.*}_corrected.txt > ${NONPAR%.*}_clean.out && mv ${NONPAR%.*}_clean.out ${CHRX%.*}_clean.out
	cat ${AUTO%.*}_clean.out <(tail -n+2 ${PAR%.*}_clean.out) > ${AUTO%.*}_tmp && mv ${AUTO%.*}_tmp ${AUTO%.*}_clean.out
fi

# Add header:
echo SNP CHR POS minor major N N_case N_control MAF MAF_case MAF_control P effect SE > ${AUTO%.*}_header.txt

# Merge the two together:
cat ${AUTO%.*}_header.txt ${AUTO%.*}_clean.out ${CHRX%.*}_clean.out > ${AUTO%_autosomes*}_wgs.out

# Remove the *_clean.out and header files:
rm ${AUTO%.*}_header.txt ${AUTO%.*}_clean.out ${AUTO%.*}_corrected.txt

if [[ -f "$CHRX" ]]; then
	rm ${CHRX%.*}_clean.out ${CHRX%.*}_corrected.txt
else
	rm ${PAR%.*}_clean.out ${PAR%.*}_corrected.txt
	rm ${CHRX%.*}_clean.out ${NONPAR%.*}_corrected.txt
fi

