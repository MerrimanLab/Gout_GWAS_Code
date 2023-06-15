#! /bin/bash

# Script to run FINEMAP for a given locus

LOCUS=${1}
SEX=${2}

sed "s/X/${LOCUS}/g" FINEMAP/finemap_master_sss_default.${SEX}.txt > FINEMAP/finemap_master_sss_default.${SEX}.${LOCUS}.txt

# Default FINEMAP, but with extra thresholds
#
# The 2.5th percentile and 97.5th percentile (i.e. 95% interval) of the effect
# size in the EUR GWAS data were -0.1134 and 0.1129, respectively. Taking the
# greater absolute value (0.1134) and using it in equation 8 in 2007 Wakefield
# paper, it gives 0.058 as the prior standard deviation of effect size with
# 0.95 probability.
./FINEMAP/finemap_v1.4_x86_64/finemap_v1.4_x86_64 \
	--sss \
	--in-files FINEMAP/finemap_master_sss_default.${SEX}.${LOCUS}.txt \
	--log \
	--flip-beta \
	--prob-cred-set 0.99 \
	--prior-std 0.058 \
	--n-conv-sss 10000 \
	--n-configs-top 1000

rm FINEMAP/finemap_master_sss_default.${SEX}.${LOCUS}.txt

