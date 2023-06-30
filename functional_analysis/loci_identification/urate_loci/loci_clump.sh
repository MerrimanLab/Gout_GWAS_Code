#! /bin/bash
# Run with parallel:
# parallel --xapply -j10 bash src/clump/loci_clump.sh dat/urate_clumping/{1}.txt {2} ::: $(cut -f7 dat/urate_clumping/crude_clumping.txt) ::: $(cut -f8 dat/urate_clumping/crude_clumping.txt)

module load plink/plink1.9b6.10

INPUT=$1
CLUMP_KB=$2
OUTPUT=$(echo ${INPUT%%.*} | sed 's/dat/res/g')

# Script to run LD-based clumping on a locus region

plink --bfile dat/1kgp/EUR_wgs.no_relatives.no_indel.biallelic --clump ${INPUT} --clump-p1 5e-8 --clump-p2 1e-5 --clump-r2 0.01 --clump-kb ${CLUMP_KB} --out ${OUTPUT} --clump-allow-overlap

