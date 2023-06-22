#! /bin/bash

# Script to list out all of the mQTL sites that are significantly altered by
# a given variant (in chr_pos).
#
# Usage:
# bash src/get_sig_meqtl.sh 4_89052323 <output_dir>
#
# Notes:
# - Ignores INDELs in the mQTL data
# - Only consider cis-mQTL

CHR=$(echo $1 | sed -e 's/_.*//g' -e 's/X/23/g')

grep -v 'INDEL' /Volumes/archive/merrimanlab/reference_files/GoDMC_meQTL/cleaned/cis/GoDMC_cis_meQTL.chr${CHR}.tsv | grep -Fw ${1} | cut -f1,8 > ${2}/${1}_mQTL.txt

