#! /bin/bash

module load bcftools/bcftools-1.11
module load htslib/htslib-1.11

bcftools view -i NSM=1 data/dbsnp_missense/GCF_000001405.25.vcf.gz -Oz -o data/dbsnp_missense/missense.vcf.gz

# Make UCSC bed format file with the missense VCF file
zcat data/dbsnp_missense/missense.vcf.gz | grep -v '^#' | cut -f1-3 | grep 'NC_' | sed -e 's/^NC_0\+/chr/g' -e 's/\.[0-9]\+//g' | awk '{print $1, $2, $2, $3} '> data/dbsnp_missense/missense.bed

