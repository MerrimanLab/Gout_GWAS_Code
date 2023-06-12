#! /bin/bash

# Merge all the log files into one

OUTPUT_DIR=$1
TRAIT=$2

cat genetic_correlation/output_header.txt > ${OUTPUT_DIR}/rg_result.tmp

for i in $(ls ${OUTPUT_DIR}/*.log); do
paste <(echo ${TRAIT} ${i##*/}) <(tail -n 4 ${i} | head -1 | tr -s ' ' '\t' | cut -f3-) | tr -s ' ' '\t' | sed 's/.log//g' >> ${OUTPUT_DIR}/rg_result.tmp
done

Rscript genetic_correlation/add_trait_info.R ${OUTPUT_DIR}/rg_result.tmp

rm ${OUTPUT_DIR}/rg_result.tmp

