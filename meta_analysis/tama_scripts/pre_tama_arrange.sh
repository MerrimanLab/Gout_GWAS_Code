#! /bin/bash

# Bash script to rearrange the merged summary stats and remove the header
# before runnint MANTRA for TAMA

OUT=${1%.*}

# Reformat the data into MANTRA input format, and substitute NAs to 0:
awk 'NR > 1 {for(i = 1;i <= NF;i++) { printf("%s%s", $i,  (i % 5 != 0) ? OFS : "\n") } }' ${1} | sed 's/NA/0/g' > ${OUT}_pre_tama.txt

