#! /bin/bash

# Script to pull out significant eQTLs from ImmuNexUT data set
# Use the same SNP list file as in OneK1K

IMMU_FILES=/Volumes/archive/merrimanlab/reference_files/ImmuNexUT/raw/*Mono*
QUERY=$1

awk 'BEGIN {OFS = "\t"}; NR == 1 {print "cell_type", $0}' /Volumes/archive/merrimanlab/reference_files/ImmuNexUT/raw/CD16p_Mono_conditional_eQTL_FDR0.05.txt > data/gene_prioritisation/immunexut_eqtl.txt
grep -Fwf ${QUERY} ${IMMU_FILES} | sed -e 's/.*\///g' -e 's/_conditional_eQTL_FDR0.05.txt:/ /g' | tr ' ' '\t' >> data/gene_prioritisation/immunexut_eqtl.txt

