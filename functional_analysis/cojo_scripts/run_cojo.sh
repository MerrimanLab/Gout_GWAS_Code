#! /bin/bash

# Script to run the whole GCTA cojo pipeline
#
# Usage:
# bash src/cojo_script/run_cojo.sh <input_file> <output_directory>
#
# Input file contain two columns:
# 1) SNP - used to label the output file
# 2) region - used to pull out the region from UKBB and GCTA cojo
#
# "region" must be in the form of "CHR:START-END" and must be padded with "0"s
# if the chromosome is a single digit (e.g. 01:10000-100000)

QUERY=$1 # SNP and region file
SUMSTATS=$2 # Summary stats that goes in GCTA cojo
OUTDIR=$3 # Output directory

if [[ -z $4 ]]; then
    CPU=5
else
    CPU=$4
fi

module load gcta/gcta_1.93.2beta

# Make a list of variants present in the EUR GWAS summary stats:
mkdir -p ${OUTDIR}
bash src/cojo_script/make_sumstat_variant_list.sh ${OUTDIR}

# Pull out bfile from UKBB:
parallel -j ${CPU} --xapply bash src/cojo_script/make_ukbb_bfile.sh {1} {2} ${OUTDIR} ::: $(cut -f2 ${QUERY}) ::: $(cut -f1 ${QUERY})

# Run GCTA cojo on each of the regions:
parallel -j ${CPU} gcta --bfile ${OUTDIR}/{1} --maf 0.01 --cojo-file ${SUMSTATS} --extract ${OUTDIR}/{1}.snplist --cojo-slct --out ${OUTDIR}/{1}.cojo_results ::: $(cut -f1 ${QUERY})

