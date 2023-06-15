#! /bin/bash

# Run SLALOM on a locus to see if locus is fit for finemapping

export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OMP_NUM_THREADS=1
export GOTO_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

SUMMARY=$1
TMP_NAME=${SUMMARY##*/}
OUTFILE=$2/${TMP_NAME%.*}

# The 2.5th percentile and 97.5th percentile (i.e. 95% interval) of the effect
# size in the EUR GWAS data were -0.1134 and 0.1129, respectively. Taking the
# greater absolute value (0.1134) and using it in equation 8 in 2007 Wakefield
# paper, it gives 0.058 as the prior standard deviation of effect size with
# 0.95 probability.

PYSPARK_SUBMIT_ARGS="--driver-memory 64g pyspark-shell" \
python src/slalom/slalom.py \
    --snp ${SUMMARY} \
    --export-r \
    --lead-variant-choice prob \
    --dentist-s \
    --abf \
    --abf-prior-variance 0.058 \
    --reference-genome GRCh37 \
    --out ${OUTFILE}.result.txt
