#! /bin/bash

# Small script to clean up the MANTRA TAMA results.

TAMA_RES=$1
CLEAN_RES=${1%.*}
RES_DIR=${1%/*}

# Add correct sample size for the variants that had more than a million
# samples:
awk ' NR == FNR {AFR[$1] = $4; EAS[$1] = $5; EUR[$1] = $6; LAT[$1] = $7; sample_size[$1] = $8; next}; FNR > 1 {$1 = $2"_"$3"_"$4"_"$5; $9 = sample_size[$1]; print $0, AFR[$1], EAS[$1], EUR[$1], LAT[$1]}  ' ${RES_DIR}/tama_sample_size.txt ${TAMA_RES} | tr -s ' ' '\t' > ${CLEAN_RES}.tmp

# Add header to the TAMA results, and keep only the variants that had
# information from all four ancestries:
# NOTE: some unique SNP IDs are truncated (likely due to fixed-width column
# output from MANTRA)
awk 'BEGIN {print "SNP", "CHR", "POS", "REF", "OTH", "N_STUD", "logBF", "PPA", "N", "DIR", "N_AFR", "N_EAS", "N_EUR", "N_LAT", "cpid"}; $6 == 4 {print $0, $2"_"$3}' ${CLEAN_RES}.tmp | tr -s ' ' '\t' > ${CLEAN_RES}.tmp2

# Clean beta:
awk 'BEGIN {print "SNP", "CHR", "POS", "N_STUD", "study1", "mean_effect1", "effect_sd1", "study2", "mean_effect2", "effect_sd2", "study3", "mean_effect3", "effect_sd3", "study4", "mean_effect4", "effect_sd4", "cpid"}; $4 == 4 {print $0, $2"_"$3}' ${TAMA_RES%%.out}.beta.out | tr -s ' ' '\t' > ${CLEAN_RES}.beta.tmp

# Generate rsIDs:
cut -f 2-3 ${CLEAN_RES}.tmp2 | sort | uniq > ${CLEAN_RES}.chrpos_info.txt
bash src/snptrack.sh ${CLEAN_RES}.chrpos_info.txt
awk 'BEGIN {print "chr_bp", "SNP"}; {print $0}' ${CLEAN_RES}.chrpos_info.result.txt > ${RES_DIR}/tmp && mv ${RES_DIR}/tmp ${CLEAN_RES}.chrpos_info.result.txt

# Merge the TAMA results, beta, and the rsID info together:
Rscript src/merge_tama.R ${CLEAN_RES}.tmp2 ${CLEAN_RES}.beta.tmp ${CLEAN_RES}.chrpos_info.result.txt && mv ${CLEAN_RES}.tmp2 ${CLEAN_RES}_clean.out

rm ${CLEAN_RES}*tmp*
