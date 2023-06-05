#! /bin/bash

# Script to remove all the N-filtered variants from the TAMA results.

FILE=$1.out
OUT=${1%.*}

# Remove variants that have been filtered out in other ancestries due to sample
# size:
#
# NOTE: no need to filter using the TAMA result itself, since any variant that
# didn't pass the single-ancestry N-filter shouldn't be included in the TAMA
# result
grep -Fwf ${FILE%.*}.npass.list ${FILE} > ${FILE%.*}.nfiltered.txt

