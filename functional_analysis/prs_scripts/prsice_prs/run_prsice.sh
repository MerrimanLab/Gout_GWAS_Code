#! /bin/bash

# Run PRSICE-2 with the EUR gout GWAS and the UKBB genotype data

prsice -b dat/prs/eur_gout_sumstats.convp.full.txt --chr CHR --A1 A1 --A2 A2 --snp SNP --bp BP --stat OR --pvalue P --target dat/prs/ukb_geno/ukbb_chr#,dat/prs/ukb_geno/ukbb_gout_pheno.fam --binary-target T --clump-kb 500 --clump-r2 0.01 --perm 10000 --out res/prs/prsice_himem -n 30 --memory 500Gb --chr-id c:L --ultra --extract res/prs/prsice_himem.valid

prsice -b dat/prs/eur_gout_sumstats.convp.male.txt --chr CHR --A1 A1 --A2 A2 --snp SNP --bp BP --stat OR --pvalue P --target dat/prs/ukb_geno/ukbb_chr#,dat/prs/ukb_geno/ukbb_gout_pheno.fam --binary-target T --clump-kb 500 --clump-r2 0.01 --perm 10000 --out res/prs/prsice_himem -n 30 --memory 500Gb --chr-id c:L --ultra --extract res/prs/prsice_himem.valid --keep dat/prs/ukbb_male_samples.txt

prsice -b dat/prs/eur_gout_sumstats.convp.female.txt --chr CHR --A1 A1 --A2 A2 --snp SNP --bp BP --stat OR --pvalue P --target dat/prs/ukb_geno/ukbb_chr#,dat/prs/ukb_geno/ukbb_gout_pheno.fam --binary-target T --clump-kb 500 --clump-r2 0.01 --perm 10000 --out res/prs/prsice_himem -n 30 --memory 500Gb --chr-id c:L --ultra --extract res/prs/prsice_himem.valid --keep dat/prs/ukbb_male_samples.txt
