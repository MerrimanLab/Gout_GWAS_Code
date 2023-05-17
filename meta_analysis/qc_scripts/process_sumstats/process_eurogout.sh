#! /bin/bash

# Script to process Eurogout data (adding MAF)

SUMMARY=$1
FREQ=$2

# Make a list of variants present in the summary stats
cut -f1 ${SUMMARY} > ${SUMMARY%.*}_tmp

# NOTE: approx. 1000 variants are not present in the frequency file, compared
# to the summary stats.
# grep -Fwf ${SUMMARY%.*}_tmp ${FREQ} | wc -l

# grep out the rsID from FREQ
grep -Fwf ${SUMMARY%.*}_tmp ${FREQ} > ${FREQ%.frq.cc.*}_reduced.freq

# Add MAF information to summary stats
Rscript src/merge_maf_plink.R ${SUMMARY} ${FREQ%.frq.cc.*}_reduced.freq

