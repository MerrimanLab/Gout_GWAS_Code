#! /bin/bash

# Script to run the smartpca program from EIGENSOFT and clean the output

POP=$1
DIR=$2

sed "s/blah/${POP}/g" ${DIR}/src/ldsc_scripts/par.smartpca > ${DIR}/src/ldsc_scripts/par.smartpca.${POP}

${DIR}/src/ldsc_scripts/EIG-6.1.4/bin/smartpca -p src/ldsc_scripts/par.smartpca.${POP} > ${DIR}/data/ldsc/1kgp_pc/${POP}_eigen.log

cat ${DIR}/data/ldsc/1kgp_pc/${POP}_eigen.pca.evec | tr -s " " "\t" | sed "s/^\t//g" | cut -f1-11 | awk 'NR > 1 {print $1, $0}' | tr "\t" " " > ${DIR}/data/ldsc/1kgp_pc/${POP}_eigen.pca.evec.clean
