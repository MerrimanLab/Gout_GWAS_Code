#! /bin/bash

# Script to run the convertf program from EIGENSOFT

POP=$1
DIR=$2

sed "s/blah/${POP}/g" ${DIR}/src/ldsc_scripts/par.PED.EIGENSTRAT > ${DIR}/src/ldsc_scripts/par.PED.EIGENSTRAT.${POP}

${DIR}/src/ldsc_scripts/EIG-6.1.4/bin/convertf -p ${DIR}/src/ldsc_scripts/par.PED.EIGENSTRAT.${POP} > ${DIR}/data/ldsc/1kgp_plink/${POP}/${POP}_convertf.log

