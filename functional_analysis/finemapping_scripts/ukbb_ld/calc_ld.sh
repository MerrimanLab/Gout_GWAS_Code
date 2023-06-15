#! /bin/bash

# Generate LD matrix for fine-mapping, given a region and rsID (for output
# naming)
# NOTE: need to make a list of variants from the summary stats to calculate LD
# with (`cut -f<rsid column> summary_stats > data/ukbb_ld/sumstats_variants.txt)

export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OMP_NUM_THREADS=1
export GOTO_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

module load bgenix/bgenix-1.1.8
module load bcftools/bcftools-1.11
module load plink/plink1.9b6.21

REGION=$(echo $1 | sed 's/^23:/X:/g')
CHR=$(echo ${REGION} | sed -e 's/:.*//g' -e 's/^0//g')
SNP=$2
UKBDIR=/Volumes/scratch/merrimanlab/ukbio/EGAD00010001474/
OUTDIR=$3

echo ${REGION} ${CHR}

# Pull out UKBBB genotypes:
if [[ $CHR == "X" ]]; then
bgenix -g ${UKBDIR}/ukb_imp_chr${CHR}_v3.bgen -i ${UKBDIR}/ukb_imp_chr${CHR}_v3.bgen.bgi -incl-range ${REGION} -vcf | bcftools reheader -h ${UKBDIR}/bgen_to_vcf/new_header_chrX.txt | bgzip -c > ${OUTDIR}/${SNP}.vcf.gz
else
bgenix -g ${UKBDIR}/ukb_imp_chr${CHR}_v3.bgen -i ${UKBDIR}/ukb_imp_chr${CHR}_v3.bgen.bgi -incl-range ${REGION} -vcf | bcftools reheader -h ${UKBDIR}/bgen_to_vcf/new_header.txt | bcftools annotate --rename-chrs ${UKBDIR}/bgen_to_vcf/rename_contigs.txt | bgzip -c > ${OUTDIR}/${SNP}.vcf.gz
fi

# Make a keep file of all the variants present in UKBB
zcat ${OUTDIR}/${SNP}.vcf.gz | grep -v '^#' | cut -f1-5 | grep -Fwf ${OUTDIR}/sumstats_variants.txt | cut -f3 > ${OUTDIR}/${SNP}.keep
zcat ${OUTDIR}/${SNP}.vcf.gz | grep -v '^#' | cut -f1-5 | grep -Ff ${OUTDIR}/sumstats_variants.cpid.txt | cut -f3 >> ${OUTDIR}/${SNP}.keep

sort ${OUTDIR}/${SNP}.keep | uniq > ${OUTDIR}/${SNP}.keep.tmp && mv ${OUTDIR}/${SNP}.keep.tmp ${OUTDIR}/${SNP}.keep

# Convert VCF to bfile, calculate MAF, and generate LD matrix
plink --vcf ${OUTDIR}/${SNP}.vcf.gz --extract ${OUTDIR}/${SNP}.keep --keep src/ukbb_ld/keep_ids.txt --maf 0.01 --make-bed --out ${OUTDIR}/${SNP}
plink --bfile ${OUTDIR}/${SNP} --pheno src/ukbb_ld/keep_ids.txt --allow-no-sex --freq --out ${OUTDIR}/${SNP}
plink --bfile ${OUTDIR}/${SNP} --r square --out ${OUTDIR}/${SNP}

# Remove VCF
rm ${OUTDIR}/${SNP}.vcf.gz

