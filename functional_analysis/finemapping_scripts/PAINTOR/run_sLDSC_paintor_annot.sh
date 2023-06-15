#! /bin/bash
export OPENBLAS_NUM_THREADS=10
export MKL_NUM_THREADS=10
export OMP_NUM_THREADS=10
export GOTO_NUM_THREADS=10
export NUMEXPR_NUM_THREADS=10

# Script to generate reference LD score for all the PAINTOR annotations in 1000
# Genomes EUR reference and run LDSC with the EUR summary stats to figure out
# which annotations to use for PAINTOR.

KGP_REF=data/ldsc/1kgp_plink/EUR/EUR_chr

# Make directory to put all the annotation files:
mkdir -p data/ldsc/1kgp_ref_paintor/ldscores/{Hnisz_Cell2013_SuperEnhancer,RoadMap_Assayed_NarrowPeak,RoadMap_Enhancers,TFBS,FANTOM5,Maurano_Science2012_DHS,Roadmap_ChromeHMM_15state,RoadMap_Imputed_NarrowPeak,Thurman_Nature2012_DHS,GeneElements_Gencode,RoadMap_Dyadic,RoadMap_Promoter}

# Convert all the bed files into tab-delimited files:
for i in {Hnisz_Cell2013_SuperEnhancer,RoadMap_Assayed_NarrowPeak,RoadMap_Enhancers,TFBS,FANTOM5,Maurano_Science2012_DHS,Roadmap_ChromeHMM_15state,RoadMap_Imputed_NarrowPeak,Thurman_Nature2012_DHS,GeneElements_Gencode,RoadMap_Dyadic,RoadMap_Promoter} ; do
	parallel -j 40 'cat {} | tr -s " " "\t" > {}.tmp && mv {}.tmp {}' ::: $(find data/ldsc/ref_files/Functional_Annotations/${i}/ | tail -n+2)
done

# Generate annotation file:
for i in {Hnisz_Cell2013_SuperEnhancer,RoadMap_Enhancers,TFBS,FANTOM5,Maurano_Science2012_DHS,Thurman_Nature2012_DHS,GeneElements_Gencode,RoadMap_Dyadic,RoadMap_Promoter} ; do
for k in {1..22} ; do
parallel -j 140 make_annot.py --bed-file {} --bimfile ${KGP_REF}${k}.mac5.bim --annot-file data/ldsc/1kgp_ref_paintor/ldscores/${i}/{/}.chr${k}.annot.gz ::: $(find data/ldsc/ref_files/Functional_Annotations/${i}/ | tail -n+2)
ls data/ldsc/1kgp_ref_paintor/ldscores/${i}/*chr1.annot.gz | sed  -e 's/^.*\///g' -e 's/.chr.*//g' | tr '\n' '\t' | sed 's/\t$/\n/g' > data/ldsc/1kgp_ref_paintor/ldscores/${i}/header.txt
gunzip data/ldsc/1kgp_ref_paintor/ldscores/${i}/*chr${k}.annot.gz
paste data/ldsc/1kgp_ref_paintor/ldscores/${i}/*chr${k}.annot | tr ' ' '\t' | tail -n+2 > data/ldsc/1kgp_ref_paintor/ldscores/${i}/${i}.chr${k}.annot.tmp
cat data/ldsc/1kgp_ref_paintor/ldscores/${i}/header.txt data/ldsc/1kgp_ref_paintor/ldscores/${i}/${i}.chr${k}.annot.tmp > data/ldsc/1kgp_ref_paintor/ldscores/${i}/${i}.chr${k}.annot && gzip data/ldsc/1kgp_ref_paintor/ldscores/${i}/${i}.chr${k}.annot && rm data/ldsc/1kgp_ref_paintor/ldscores/${i}/*chr${k}.annot data/ldsc/1kgp_ref_paintor/ldscores/${i}/*chr${k}.annot.tmp
done
done

# Separate loop for annotations that have too many annotations to paste together
for i in {RoadMap_Assayed_NarrowPeak,Roadmap_ChromeHMM_15state,RoadMap_Imputed_NarrowPeak} ; do
for k in {1..22} ; do
parallel -j 200 make_annot.py --bed-file {} --bimfile ${KGP_REF}${k}.mac5.bim --annot-file data/ldsc/1kgp_ref_paintor/ldscores/${i}/{/}.chr${k}.annot.gz ::: $(find data/ldsc/ref_files/Functional_Annotations/${i}/ | tail -n+2)
ls data/ldsc/1kgp_ref_paintor/ldscores/${i}/*chr${k}.annot.gz | sed  -e 's/^.*\///g' -e 's/.chr.*//g' | tr '\n' '\t' | sed 's/\t$/\n/g' > data/ldsc/1kgp_ref_paintor/ldscores/${i}/header.txt
gunzip data/ldsc/1kgp_ref_paintor/ldscores/${i}/*chr${k}.annot.gz
parallel "paste data/ldsc/1kgp_ref_paintor/ldscores/${i}/E{}[_-]*chr${k}.annot | tr ' ' '\t' | tail -n+2 > data/ldsc/1kgp_ref_paintor/ldscores/${i}/${i}.E{}.chr${k}.annot.tmp " ::: {001..129}
paste data/ldsc/1kgp_ref_paintor/ldscores/${i}/${i}.E*.chr${k}.annot.tmp > data/ldsc/1kgp_ref_paintor/ldscores/${i}/${i}.chr${k}.annot.tmp
cat data/ldsc/1kgp_ref_paintor/ldscores/${i}/header.txt data/ldsc/1kgp_ref_paintor/ldscores/${i}/${i}.chr${k}.annot.tmp > data/ldsc/1kgp_ref_paintor/ldscores/${i}/${i}.chr${k}.annot && gzip data/ldsc/1kgp_ref_paintor/ldscores/${i}/${i}.chr${k}.annot && rm data/ldsc/1kgp_ref_paintor/ldscores/${i}/*chr${k}.annot data/ldsc/1kgp_ref_paintor/ldscores/${i}/*chr${k}.annot.tmp
done
done

# Paste together the original bim file with the annot file to make a proper
# annotation file (not thin annotation file)
for i in {1..22} ; do
	awk 'BEGIN {print "CHR", "BP", "SNP", "CM", "base"}; {print $1, $4, $2, $3, 1}' data/ldsc/1kgp_plink/EUR/EUR_chr${i}.mac5.bim | tr ' ' '\t' > data/ldsc/1kgp_ref_paintor/chr${i}.bim
done

parallel "gunzip data/ldsc/1kgp_ref_paintor/ldscores/{2}/{2}.chr{1}.annot.gz && paste data/ldsc/1kgp_ref_paintor/chr{1}.bim data/ldsc/1kgp_ref_paintor/ldscores/{2}/{2}.chr{1}.annot > data/ldsc/1kgp_ref_paintor/ldscores/{2}/{2}.chr{1}.annot.tmp && mv data/ldsc/1kgp_ref_paintor/ldscores/{2}/{2}.chr{1}.annot.tmp data/ldsc/1kgp_ref_paintor/ldscores/{2}/{2}.chr{1}.annot && gzip data/ldsc/1kgp_ref_paintor/ldscores/{2}/{2}.chr{1}.annot" ::: {1..22} ::: {Hnisz_Cell2013_SuperEnhancer,RoadMap_Assayed_NarrowPeak,RoadMap_Enhancers,TFBS,FANTOM5,Maurano_Science2012_DHS,Roadmap_ChromeHMM_15state,RoadMap_Imputed_NarrowPeak,Thurman_Nature2012_DHS,GeneElements_Gencode,RoadMap_Dyadic,RoadMap_Promoter}

# Make reference LDSC file with the annotations
parallel -j 80 cov-ldsc.py --bfile data/ldsc/1kgp_plink/EUR/EUR_chr{1}.mac5 --annot data/ldsc/1kgp_ref_paintor/ldscores/{2}/{2}.chr{1}.annot.gz --cov data/ldsc/1kgp_pc/EUR_eigen.pca.evec.clean --print-snps data/ldsc/ref_files/w_hm3_nomhc_snplist_plink.txt --l2 --ld-wind-cm 20 --out data/ldsc/1kgp_ref_paintor/ldscores/{2}/{2}.chr{1} ::: {1..22} ::: {Hnisz_Cell2013_SuperEnhancer,RoadMap_Enhancers,TFBS,FANTOM5,Maurano_Science2012_DHS,Thurman_Nature2012_DHS,GeneElements_Gencode,RoadMap_Dyadic,RoadMap_Promoter,RoadMap_Assayed_NarrowPeak,RoadMap_Imputed_NarrowPeak,Roadmap_ChromeHMM_15state}

# Generate ldcts files for all the annotations

find data/ldsc/1kgp_ref_paintor/ldscores/TFBS/ -name "*chr3.annot.gz" | sed 's/chr3.*//g' | awk '{print $1, $1"chr"}' > data/ldsc/1kgp_ref_paintor/ldscores/TFBS/TFBS.ldcts

find data/ldsc/1kgp_ref_paintor/ldscores/Hnisz_Cell2013_SuperEnhancer/ -name "*chr3.annot.gz" | sed 's/chr3.*//g' | awk '{print $1, $1"chr"}' > data/ldsc/1kgp_ref_paintor/ldscores/Hnisz_Cell2013_SuperEnhancer/Hnisz_Cell2013_SuperEnhancer.ldcts

find data/ldsc/1kgp_ref_paintor/ldscores/RoadMap_Enhancers/ -name "*chr3.annot.gz" | sed 's/chr3.*//g' | awk '{print $1, $1"chr"}' > data/ldsc/1kgp_ref_paintor/ldscores/RoadMap_Enhancers/RoadMap_Enhancers.ldcts

find data/ldsc/1kgp_ref_paintor/ldscores/FANTOM5/ -name "*chr3.annot.gz" | sed 's/chr3.*//g' | awk '{print $1, $1"chr"}' > data/ldsc/1kgp_ref_paintor/ldscores/FANTOM5/FANTOM5.ldcts

find data/ldsc/1kgp_ref_paintor/ldscores/RoadMap_Dyadic/ -name "*chr3.annot.gz" | sed 's/chr3.*//g' | awk '{print $1, $1"chr"}' > data/ldsc/1kgp_ref_paintor/ldscores/RoadMap_Dyadic/RoadMap_Dyadic.ldcts

find data/ldsc/1kgp_ref_paintor/ldscores/RoadMap_Promoter/ -name "*chr3.annot.gz" | sed 's/chr3.*//g' | awk '{print $1, $1"chr"}' > data/ldsc/1kgp_ref_paintor/ldscores/RoadMap_Promoter/RoadMap_Promoter.ldcts

find data/ldsc/1kgp_ref_paintor/ldscores/Roadmap_ChromeHMM_15state/ -name "*chr3.annot.gz" | sed 's/chr3.*//g' | awk '{print $1, $1"chr"}' > data/ldsc/1kgp_ref_paintor/ldscores/Roadmap_ChromeHMM_15state/Roadmap_ChromeHMM_15state.ldcts

# For Hnisz annotations, need to remove BI_CD8_Naive_7pool.temp needs to be removed, as it causes singular matrix error
find data/ldsc/1kgp_ref_paintor/ldscores/Hnisz_Cell2013_SuperEnhancer/ -name "*chr3.annot.gz" | sed 's/chr3.*//g' | awk '{print $1, $1"chr"}' | grep -v 'temp' > data/ldsc/1kgp_ref_paintor/ldscores/Hnisz_Cell2013_SuperEnhancer/Hnisz_Cell2013_SuperEnhancer.ldcts

# all_fdr0.05_hot annotation causes singular matrix error
find data/ldsc/1kgp_ref_paintor/ldscores/Maurano_Science2012_DHS/ -name "*chr3.annot.gz" | sed 's/chr3.*//g' | awk '{print $1, $1"chr"}' | grep -v 'all_fdr' > data/ldsc/1kgp_ref_paintor/ldscores/Maurano_Science2012_DHS/Maurano_Science2012_DHS.ldcts

find data/ldsc/1kgp_ref_paintor/ldscores/Thurman_Nature2012_DHS/ -name "*chr3.annot.gz" | sed 's/chr3.*//g' | awk '{print $1, $1"chr"}' > data/ldsc/1kgp_ref_paintor/ldscores/Thurman_Nature2012_DHS/Thurman_Nature2012_DHS.ldcts

# Selenocysteine annotation causes singular matrix error (only 0s in the annotation)
find data/ldsc/1kgp_ref_paintor/ldscores/GeneElements_Gencode/ -name "*chr3.annot.gz" | sed 's/chr3.*//g' | awk '{print $1, $1"chr"}' | grep -v 'Selenocysteine' > data/ldsc/1kgp_ref_paintor/ldscores/GeneElements_Gencode/GeneElements_Gencode.ldcts

# Somehow there are a couple of temporary files in the list, so remove them.
# There was a full annotation file (RoadMap_Assayed_NarrowPeak.ldcts) included
# in the ldcts file - remove it afterwards
find data/ldsc/1kgp_ref_paintor/ldscores/RoadMap_Assayed_NarrowPeak/ -name "*chr3.annot.gz" | sed 's/chr3.*//g' | awk '{print $1, $1"chr"}' | grep -v '\/\.E' > data/ldsc/1kgp_ref_paintor/ldscores/RoadMap_Assayed_NarrowPeak/RoadMap_Assayed_NarrowPeak.ldcts

# These tissues/histone marks caused linalg error:
find data/ldsc/1kgp_ref_paintor/ldscores/RoadMap_Imputed_NarrowPeak/ -name "*chr3.annot.gz" | sed 's/chr3.*//g' | awk '{print $1, $1"chr"}' | grep -vf data/ldsc/1kgp_ref_paintor/ldscores/RoadMap_Imputed_NarrowPeak/remove_linalg_annot.txt > data/ldsc/1kgp_ref_paintor/ldscores/RoadMap_Imputed_NarrowPeak/RoadMap_Imputed_NarrowPeak.ldcts

# Run LD score regression using the .ldcts files

mkdir -p results/ldsc/paintor_ldsc

# Use baseline model to figure out relevant annotations
cov-ldsc.py --h2-cts data/fema_data/EUR/EUR_meta_full1_clean_rsid.nfiltered.sumstats.gz --ref-ld-chr data/ldsc/1kgp_ref_ldsc/EUR/EUR_chr --ref-ld-chr-cts data/ldsc/1kgp_ref_paintor/ldscores/TFBS/TFBS.ldcts --w-ld-chr data/ldsc/1kgp_ref_ldsc/EUR/EUR_weights_chr --out results/ldsc/paintor_ldsc/TFBS.baseline

cov-ldsc.py --h2-cts data/fema_data/EUR/EUR_meta_full1_clean_rsid.nfiltered.sumstats.gz --ref-ld-chr data/ldsc/1kgp_ref_ldsc/EUR/EUR_chr --ref-ld-chr-cts data/ldsc/1kgp_ref_paintor/ldscores/Hnisz_Cell2013_SuperEnhancer/Hnisz_Cell2013_SuperEnhancer.ldcts --w-ld-chr data/ldsc/1kgp_ref_ldsc/EUR/EUR_weights_chr --out results/ldsc/paintor_ldsc/Hnisz_Cell2013_SuperEnhancer.baseline

cov-ldsc.py --h2-cts data/fema_data/EUR/EUR_meta_full1_clean_rsid.nfiltered.sumstats.gz --ref-ld-chr data/ldsc/1kgp_ref_ldsc/EUR/EUR_chr --ref-ld-chr-cts data/ldsc/1kgp_ref_paintor/ldscores/RoadMap_Enhancers/RoadMap_Enhancers.ldcts --w-ld-chr data/ldsc/1kgp_ref_ldsc/EUR/EUR_weights_chr --out results/ldsc/paintor_ldsc/RoadMap_Enhancers.baseline

cov-ldsc.py --h2-cts data/fema_data/EUR/EUR_meta_full1_clean_rsid.nfiltered.sumstats.gz --ref-ld-chr data/ldsc/1kgp_ref_ldsc/EUR/EUR_chr --ref-ld-chr-cts data/ldsc/1kgp_ref_paintor/ldscores/FANTOM5/FANTOM5.ldcts --w-ld-chr data/ldsc/1kgp_ref_ldsc/EUR/EUR_weights_chr --out results/ldsc/paintor_ldsc/FANTOM5.baseline

cov-ldsc.py --h2-cts data/fema_data/EUR/EUR_meta_full1_clean_rsid.nfiltered.sumstats.gz --ref-ld-chr data/ldsc/1kgp_ref_ldsc/EUR/EUR_chr --ref-ld-chr-cts data/ldsc/1kgp_ref_paintor/ldscores/Maurano_Science2012_DHS/Maurano_Science2012_DHS.ldcts --w-ld-chr data/ldsc/1kgp_ref_ldsc/EUR/EUR_weights_chr --out results/ldsc/paintor_ldsc/Maurano_Science2012_DHS.baseline

cov-ldsc.py --h2-cts data/fema_data/EUR/EUR_meta_full1_clean_rsid.nfiltered.sumstats.gz --ref-ld-chr data/ldsc/1kgp_ref_ldsc/EUR/EUR_chr --ref-ld-chr-cts data/ldsc/1kgp_ref_paintor/ldscores/Thurman_Nature2012_DHS/Thurman_Nature2012_DHS.ldcts --w-ld-chr data/ldsc/1kgp_ref_ldsc/EUR/EUR_weights_chr --out results/ldsc/paintor_ldsc/Thurman_Nature2012_DHS.baseline

cov-ldsc.py --h2-cts data/fema_data/EUR/EUR_meta_full1_clean_rsid.nfiltered.sumstats.gz --ref-ld-chr data/ldsc/1kgp_ref_ldsc/EUR/EUR_chr --ref-ld-chr-cts data/ldsc/1kgp_ref_paintor/ldscores/GeneElements_Gencode/GeneElements_Gencode.ldcts --w-ld-chr data/ldsc/1kgp_ref_ldsc/EUR/EUR_weights_chr --out results/ldsc/paintor_ldsc/GeneElements_Gencode.baseline

cov-ldsc.py --h2-cts data/fema_data/EUR/EUR_meta_full1_clean_rsid.nfiltered.sumstats.gz --ref-ld-chr data/ldsc/1kgp_ref_ldsc/EUR/EUR_chr --ref-ld-chr-cts data/ldsc/1kgp_ref_paintor/ldscores/RoadMap_Dyadic/RoadMap_Dyadic.ldcts --w-ld-chr data/ldsc/1kgp_ref_ldsc/EUR/EUR_weights_chr --out results/ldsc/paintor_ldsc/RoadMap_Dyadic.baseline

cov-ldsc.py --h2-cts data/fema_data/EUR/EUR_meta_full1_clean_rsid.nfiltered.sumstats.gz --ref-ld-chr data/ldsc/1kgp_ref_ldsc/EUR/EUR_chr --ref-ld-chr-cts data/ldsc/1kgp_ref_paintor/ldscores/RoadMap_Promoter/RoadMap_Promoter.ldcts --w-ld-chr data/ldsc/1kgp_ref_ldsc/EUR/EUR_weights_chr --out results/ldsc/paintor_ldsc/RoadMap_Promoter.baseline

cov-ldsc.py --h2-cts data/fema_data/EUR/EUR_meta_full1_clean_rsid.nfiltered.sumstats.gz --ref-ld-chr data/ldsc/1kgp_ref_ldsc/EUR/EUR_chr --ref-ld-chr-cts data/ldsc/1kgp_ref_paintor/ldscores/RoadMap_Assayed_NarrowPeak/RoadMap_Assayed_NarrowPeak.ldcts --w-ld-chr data/ldsc/1kgp_ref_ldsc/EUR/EUR_weights_chr --out results/ldsc/paintor_ldsc/RoadMap_Assayed_NarrowPeak.baseline

cov-ldsc.py --h2-cts data/fema_data/EUR/EUR_meta_full1_clean_rsid.nfiltered.sumstats.gz --ref-ld-chr data/ldsc/1kgp_ref_ldsc/EUR/EUR_chr --ref-ld-chr-cts data/ldsc/1kgp_ref_paintor/ldscores/Roadmap_ChromeHMM_15state/Roadmap_ChromeHMM_15state.ldcts --w-ld-chr data/ldsc/1kgp_ref_ldsc/EUR/EUR_weights_chr --out results/ldsc/paintor_ldsc/Roadmap_ChromeHMM_15state.baseline

cov-ldsc.py --h2-cts data/fema_data/EUR/EUR_meta_full1_clean_rsid.nfiltered.sumstats.gz --ref-ld-chr data/ldsc/1kgp_ref_ldsc/EUR/EUR_chr --ref-ld-chr-cts data/ldsc/1kgp_ref_paintor/ldscores/RoadMap_Imputed_NarrowPeak/RoadMap_Imputed_NarrowPeak.ldcts --w-ld-chr data/ldsc/1kgp_ref_ldsc/EUR/EUR_weights_chr --out results/ldsc/paintor_ldsc/RoadMap_Imputed_NarrowPeak.baseline

