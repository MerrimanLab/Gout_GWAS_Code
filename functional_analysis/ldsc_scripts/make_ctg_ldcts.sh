#! /bin/bash

# Script to generate a .ldcts file for cell-type group and cell-type specific
# LDSC analysis.

NAMES=data/ldsc/ld_annot/ct_and_ctg_bedfiles/cell_type_group_specific/names.txt

for ancestry in {EUR,EAS,AFR,LAT} ; do
	for k in {1..10} ; do
		echo data/ldsc/1kgp_ref_ctg/${ancestry}/ctg/${ancestry}_ctg. >> data/ldsc/1kgp_ref_ctg/${ancestry}/tmp
	done
	paste -d '' data/ldsc/1kgp_ref_ctg/${ancestry}/tmp <(tail -n+2 ${NAMES} | sed 's/.bed//g') | awk '{print $2, $1"."}' | tr ' ' '\t' > data/ldsc/1kgp_ref_ctg/${ancestry}/${ancestry}_ctg.ldcts
	rm data/ldsc/1kgp_ref_ctg/${ancestry}/tmp
done

