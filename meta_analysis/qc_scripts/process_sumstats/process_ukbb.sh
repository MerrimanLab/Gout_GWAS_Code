#! /bin/bash

# Script to process UKBB data

SUMMARY=$1

cut -d',' -f1-3,6-13,19-20,23- ${SUMMARY} | tr ',' ' ' > ${SUMMARY%.*}.tmp

# Double-check HWE and INFO filters:
awk '$12 > 1e-6 || $13 >= 0.3 {$6 = log($6); $8 = (10 ^ (-1 * $8)); print $0}' ${SUMMARY%.*}.tmp > ${SUMMARY%.*}_clean.txt

rm ${SUMMARY%.*}.tmp
