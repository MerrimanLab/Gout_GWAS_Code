#! /bin/bash
# Short script to clean up the Korean summary data

KOR_DAT=$1
OUT_PATH=${1%/*}

# First remove remove unwanted space from the header:
sed 's/KARE\ /KARE_/g' ${KOR_DAT} | tr ' ' '\t' > ${OUT_PATH}/korean_clean_header.txt

# Take a look at the columns:
# head -1 ${OUT_PATH}/korean_clean_header.txt | tr '\t' '\n' | nl

# Cut out the relevant columns from the summary stats, clean the header, and
# convert X to 23 and PAR to 25:
cut -f1-3,7-15 ${OUT_PATH}/korean_clean_header.txt | sed -e 's/\.b37//g' -e 's/KARE_8842://g' -e 's/pheno\.AS1_PdGt://g' -e 's/^X/23/g' -e 's/^PAR/25/g' > ${OUT_PATH}/tmp

# Log the OR and append the column:
awk 'NR == 1 {print $0, "effect", "N", "N_case", "N_control"}; NR > 1 {print $0, log($11), 8834, 520, 8314}' ${OUT_PATH}/tmp | tr ' ' '\t' > ${OUT_PATH}/korean_gout.txt

rm ${OUT_PATH}/tmp ${OUT_PATH}/korean_clean_header.txt
