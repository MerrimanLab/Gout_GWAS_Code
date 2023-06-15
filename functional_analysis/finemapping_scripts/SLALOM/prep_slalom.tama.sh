#! /bin/bash

export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OMP_NUM_THREADS=1
export GOTO_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

SEX=$1

# Pull out TAMA loci that aren't present in EUR
# (Note that even though I'm only pulling out only the TAMA region, SLALOM is
# still run on all loci)

Rscript src/slalom/pull_tama_loci.R ${SEX}

# SLALOM assumes allele2 as OTH/effect/minor allele and allele1 as REF allele.
# Need to make sure the REF allele is actually REF allele (done by checking the
# Neale's UKBB variant list)

parallel -j10 "Rscript src/slalom/check_ref_allele.R {}" ::: $(ls data/slalom_data/${SEX}/TAMA/*.clean.txt)

