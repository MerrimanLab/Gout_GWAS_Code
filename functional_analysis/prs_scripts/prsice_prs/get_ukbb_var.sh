#! /bin/bash

# Mini-script to be run with src/prs/parallel_ukbb_var.sh to pull out specific
# variants from UKBB

module load plink/plink1.9b6.10
module load bcftools/bcftools-1.11
module load bgenix/bgenix-1.1.8

UKBDIR=/Volumes/scratch/merrimanlab/ukbio/EGAD00010001474/
CHR=$1
VAR_FILE=$2

bgenix -g ${UKBDIR}/ukb_imp_chr${CHR}_v3.bgen -i ${UKBDIR}/ukb_imp_chr${CHR}_v3.bgen.bgi -incl-range ${VAR_FILE} -vcf | bcftools reheader -h ${UKBDIR}/bgen_to_vcf/new_header.txt | bgzip -c > dat/prs/ukb_geno/${VAR_FILE##*/}.vcf.gz

plink --vcf dat/prs/ukb_geno/${VAR_FILE##*/}.vcf.gz --keep src/prs/keep_ids.txt --maf 0.01 --hwe 0.000001 --snps-only just-acgt --make-bed --out dat/prs/ukb_geno/${VAR_FILE##*/}

rm dat/prs/ukb_geno/${VAR_FILE##*/}.vcf.gz

