#! /bin/bash

# Script to generate a .ldcts file for cell-type group and cell-type specific
# LDSC analysis.

NAMES=data/ldsc/ld_annot/ct_and_ctg_bedfiles/cell_type_specific/names.txt

for ancestry in {EUR,EAS,AFR,LAT} ; do
	for k in {1..220} ; do
		echo data/ldsc/1kgp_ref_ctg/${ancestry}/cts/${ancestry}_cts. >> data/ldsc/1kgp_ref_ctg/${ancestry}/tmp
	done
	paste -d '' data/ldsc/1kgp_ref_ctg/${ancestry}/tmp <(tail -n+2 ${NAMES} | cut -f2,4-5 | sed 's/_peaks.bed//g' | awk '{print $2".", $1"_"$3}' ) | awk '{print $2, $1}' | tr ' ' '\t' > data/ldsc/1kgp_ref_ctg/${ancestry}/${ancestry}_cts.ldcts
	rm data/ldsc/1kgp_ref_ctg/${ancestry}/tmp
done

