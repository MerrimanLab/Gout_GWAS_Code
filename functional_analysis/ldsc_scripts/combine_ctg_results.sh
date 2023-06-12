#! /bin/bash

# Script to combine the results from cell-type group partitioned heritability
# analysis

ANCESTRY=$1
CTS_FILE=data/ldsc/1kgp_ref_ctg/${ANCESTRY}/${ANCESTRY}_ctg.ldcts
PREFIX=$2

for i in {1..10} ; do
	head -2 ${PREFIX}.ctg${i}.results | tail -n1 | cut -f2- >> ${PREFIX}.tmp
done

paste <(cut -f1 ${CTS_FILE}) ${PREFIX}.tmp > ${PREFIX}.tmp2

cat <(head -1 ${PREFIX}.ctg1.results) ${PREFIX}.tmp2 > ${PREFIX}.cell_type_results.txt

rm ${PREFIX}.tmp*

