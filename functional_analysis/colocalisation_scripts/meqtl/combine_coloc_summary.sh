#! /bin/bash

# Just cat together all the summary results from coloc
#
# Usage:
# bash src/combine_coloc_summary.sh <path_to_results_directory>

cat <(head -1 -q ${1}/*coloc_summary.txt | sort | uniq) <(tail -q -n+2 ${1}/*coloc_summary.txt) | tr ' ' '\t' > ${1}/coloc_results.txt

