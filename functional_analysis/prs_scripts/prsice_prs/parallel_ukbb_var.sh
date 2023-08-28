#! /bin/bash

# Script to pull out a list of (somewhat) filtered UKBB variants for PRSice-2

module load plink/plink1.9b6.10

VAR_FILE=dat/prs/chr${1}_varlist.txt
CHR=$(echo ${1} | sed 's/^0//g')

split -l 4999 -d -a 3 ${VAR_FILE} dat/prs/x

parallel -j 10 bash src/prs/get_ukbb_var.sh ${CHR} {} ::: $(ls dat/prs/x*)

ls dat/prs/ukb_geno/x*bim | sed 's/.bim//g' > dat/prs/ukb_geno/chr${CHR}_mergelist.txt

plink --merge-list dat/prs/ukb_geno/chr${CHR}_mergelist.txt --make-bed --out dat/prs/ukb_geno/ukbb_chr${CHR}

rm dat/prs/ukb_geno/x* dat/prs/x* dat/prs/ukb_geno/chr${CHR}_mergelist.txt

