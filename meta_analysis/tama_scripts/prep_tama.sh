#!/bin/bash
# Script to calculate distance between populations included in TAMA
# NOTE: you need `dmatcal` binary in your PATH

COHORT=$1
DIR=data/tama_dat/${COHORT}/
INFILE=${COHORT}_merged_dat_pre_tama.txt
PREFILE=${COHORT}_merged_dat.tsv

cd ${DIR}

head -1 ${PREFILE} | tr -s '[:blank:]' '\n' | nl | grep 'presence' | sed 's/.*presence\.//g' > mantra.in

mv ${INFILE} mantra.dat

dmatcal

tar -czvf ${COHORT}_tama_dat.tar.gz mantra.dat dmat.out mantra.in

mv mantra.dat ${INFILE}

cd -

