#! /bin/bash

# Pull out regions from UKBB and process it for GCTA cojo

module load bgenix/bgenix-1.1.8
module load bcftools/bcftools-1.11
module load plink/plink1.9b6.21

REGION=$(echo $1 | sed 's/^23:/X:/g')
CHR=$(echo ${REGION} | sed -e 's/:.*//g' -e 's/^0//g')
SNP=$2
UKBDIR=/Volumes/scratch/merrimanlab/ukbio/EGAD00010001474/
OUTDIR=$3

# Make output directory, if not already there:
mkdir -p ${OUTDIR}

# Pull out UKBBB genotypes:
if [[ $CHR == "X" ]]; then
bgenix -g ${UKBDIR}/ukb_imp_chr${CHR}_v3.bgen -i ${UKBDIR}/ukb_imp_chr${CHR}_v3.bgen.bgi -incl-range ${REGION} -vcf | bcftools reheader -h ${UKBDIR}/bgen_to_vcf/new_header_chrX.txt | bgzip -c > ${OUTDIR}/${SNP}.vcf.gz
else
bgenix -g ${UKBDIR}/ukb_imp_chr${CHR}_v3.bgen -i ${UKBDIR}/ukb_imp_chr${CHR}_v3.bgen.bgi -incl-range ${REGION} -vcf | bcftools reheader -h ${UKBDIR}/bgen_to_vcf/new_header.txt | bcftools annotate --rename-chrs ${UKBDIR}/bgen_to_vcf/rename_contigs.txt | bgzip -c > ${OUTDIR}/${SNP}.vcf.gz
fi

# Make a keep file of all the variants (in terms of chr/pos) present in UKBB
zcat ${OUTDIR}/${SNP}.vcf.gz | grep -v '^#' | cut -f1-5 | grep -Fwf ${OUTDIR}/sumstat_variant_list.txt | cut -f3 > ${OUTDIR}/${SNP}.keep
zcat ${OUTDIR}/${SNP}.vcf.gz | grep -v '^#' | cut -f1-5 | grep -Ff ${OUTDIR}/sumstat_variant_list.cpid.txt | cut -f3 >> ${OUTDIR}/${SNP}.keep

sort ${OUTDIR}/${SNP}.keep | uniq > ${OUTDIR}/${SNP}.keep.tmp && mv ${OUTDIR}/${SNP}.keep.tmp ${OUTDIR}/${SNP}.keep

# Convert VCF to bfile:
plink --vcf ${OUTDIR}/${SNP}.vcf.gz --extract ${OUTDIR}/${SNP}.keep --keep /Volumes/scratch/merrimanlab/murray/ukbb_ld/keep_ids.txt --make-bed --out ${OUTDIR}/${SNP}

# Remove VCF
rm ${OUTDIR}/${SNP}.vcf.gz

# First generate allelic value IDs for the bim file, then rename the bim file
# with this new ID
Rscript src/cojo_script/make_avid.bim.R ${OUTDIR}/${SNP}.bim ${OUTDIR}/${SNP}.avid.txt

plink --bfile ${OUTDIR}/${SNP} --update-name ${OUTDIR}/${SNP}.avid.txt --make-bed --out ${OUTDIR}/${SNP}

# Now, pull out UKBB variants based on the new ID (this should remove any
# multi-allelic variants)
grep -Fwf ${OUTDIR}/sumstats_avid_only.txt ${OUTDIR}/${SNP}.bim | tr  -s ' ' '\t' | cut -f2 > ${OUTDIR}/${SNP}.avid.keep
plink --bfile ${OUTDIR}/${SNP} --extract ${OUTDIR}/${SNP}.avid.keep --make-bed --out ${OUTDIR}/${SNP}

# Using the summary stats allelic value ID file, rename the bfiles again so
# the IDs are rsIDs from the summary stats
plink --bfile ${OUTDIR}/${SNP} --update-name ${OUTDIR}/sumstats_avid.txt 1 2 --make-bed --out ${OUTDIR}/${SNP}

cut -f2 ${OUTDIR}/${SNP}.bim > ${OUTDIR}/${SNP}.snplist

rm ${OUTDIR}/${SNP}.*~
