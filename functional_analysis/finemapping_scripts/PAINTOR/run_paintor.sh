#! /bin/bash
export OPENBLAS_NUM_THREADS=10
export MKL_NUM_THREADS=10
export OMP_NUM_THREADS=10
export GOTO_NUM_THREADS=10
export NUMEXPR_NUM_THREADS=10

module load PAINTOR/PAINTOR_V3.0

mkdir -p results/credible_set/PAINTOR/output/{full,male,female}
mkdir -p results/credible_set/PAINTOR/output_multi/{full,male,female}

# Run PAINTOR with top 5 annotations

# First make a file with paths to the bed files for the significant annotations
sed -e 's:data/ldsc/1kgp_ref_paintor/tmp_annot/:data/ldsc/ref_files/Functional_Annotations/:g' -e 's/.$//g' results/ldsc/paintor_ldsc/sig_annot.uncorrelated.txt > results/ldsc/paintor_ldsc/annot_bed_paths.txt
echo data/dbsnp_missense/missense.bed >> results/ldsc/paintor_ldsc/annot_bed_paths.txt

module load PAINTOR/PAINTOR_V3.0

# Generate annotation file for each locus
parallel -j15 python2.7 src/PAINTOR_V3.0/PAINTOR_Utilities/AnnotateLocus.py --input results/ldsc/paintor_ldsc/annot_bed_paths.txt --locus data/fine_mapping/paintor_dat/full/{}.for_annot.txt --out data/fine_mapping/paintor_dat/full/{}.annot.txt --chr CHR --pos BP ::: $(cut -f1 data/ukbb_ld/full/loci_regions.txt)
parallel -j15 python2.7 src/PAINTOR_V3.0/PAINTOR_Utilities/AnnotateLocus.py --input results/ldsc/paintor_ldsc/annot_bed_paths.txt --locus data/fine_mapping/paintor_dat/male/{}.for_annot.txt --out data/fine_mapping/paintor_dat/male/{}.annot.txt --chr CHR --pos BP ::: $(cut -f1 data/ukbb_ld/male/loci_regions.txt)
parallel -j15 python2.7 src/PAINTOR_V3.0/PAINTOR_Utilities/AnnotateLocus.py --input results/ldsc/paintor_ldsc/annot_bed_paths.txt --locus data/fine_mapping/paintor_dat/female/{}.for_annot.txt --out data/fine_mapping/paintor_dat/female/{}.annot.txt --chr CHR --pos BP ::: $(cut -f1 data/ukbb_ld/female/loci_regions.txt)

# Rename the header so it's easier to refer to them in the command:
for i in $(ls data/fine_mapping/paintor_dat/*/*[0-9].annot.txt); do
cat <(echo A{1..5}) <(tail -n+2 ${i}) > ${i%%/*}/tmp && mv ${i%%/*}/tmp ${i}
done

# Rename the summary stats and annotations files:
for i in $(ls data/fine_mapping/paintor_dat/*/*.clean.txt); do
cp ${i} ${i%%.clean.txt}
done

for i in $(ls data/fine_mapping/paintor_dat/*/*.annot.txt); do
mv ${i} ${i%%.annot.txt}.annotations
done

# Run PAINTOR

parallel -j20 "echo {} > data/fine_mapping/paintor_dat/female/{}.in && PAINTOR -input data/fine_mapping/paintor_dat/female/{}.in -in data/fine_mapping/paintor_dat/female -annotations A1,A2,A3,A4,A5 -Zhead z -LDname ld -Gname {}.est.enrichment -Lname {}.logBF -enumerate 1 -out results/credible_set/PAINTOR/output/female > results/credible_set/PAINTOR/output/female/{}.paintor.log" ::: $(cut -f1 data/ukbb_ld/female/loci_regions.txt)
parallel -j80 "echo {} > data/fine_mapping/paintor_dat/male/{}.in && PAINTOR -input data/fine_mapping/paintor_dat/male/{}.in -in data/fine_mapping/paintor_dat/male -annotations A1,A2,A3,A4,A5 -Zhead z -LDname ld -Gname {}.est.enrichment -Lname {}.logBF -enumerate 1 -out results/credible_set/PAINTOR/output/male > results/credible_set/PAINTOR/output/male/{}.paintor.log" ::: $(cut -f1 data/ukbb_ld/male/loci_regions.txt)
parallel -j80 "echo {} > data/fine_mapping/paintor_dat/full/{}.in && PAINTOR -input data/fine_mapping/paintor_dat/full/{}.in -in data/fine_mapping/paintor_dat/full -annotations A1,A2,A3,A4,A5 -Zhead z -LDname ld -Gname {}.est.enrichment -Lname {}.logBF -enumerate 1 -out results/credible_set/PAINTOR/output/full > results/credible_set/PAINTOR/output/full/{}.paintor.log" ::: $(cut -f1 data/ukbb_ld/full/loci_regions.txt)

# Run specific loci with different max_causal option when the loci might have
# multiple signals:

parallel "PAINTOR -input data/fine_mapping/paintor_dat/female/{}.in -in data/fine_mapping/paintor_dat/female -annotations A1,A2,A3,A4,A5 -Zhead z -LDname ld -Gname {}.est.enrichment -Lname {}.logBF -enumerate 3 -out results/credible_set/PAINTOR/output_multi/female > results/credible_set/PAINTOR/output_multi/female/{}.paintor.log" ::: $(cut -f4 results/credible_set/PAINTOR/output_multi/female/multi_signal_loci.txt)

parallel "PAINTOR -input data/fine_mapping/paintor_dat/full/{}.in -in data/fine_mapping/paintor_dat/full -annotations A1,A2,A3,A4,A5 -Zhead z -LDname ld -Gname {}.est.enrichment -Lname {}.logBF -enumerate 3 -out results/credible_set/PAINTOR/output_multi/full > results/credible_set/PAINTOR/output_multi/full/{}.paintor.log" ::: $(cut -f4 results/credible_set/PAINTOR/output_multi/full/multi_signal_loci.txt)

parallel "PAINTOR -input data/fine_mapping/paintor_dat/male/{}.in -in data/fine_mapping/paintor_dat/male -annotations A1,A2,A3,A4,A5 -Zhead z -LDname ld -Gname {}.est.enrichment -Lname {}.logBF -enumerate 3 -out results/credible_set/PAINTOR/output_multi/male > results/credible_set/PAINTOR/output_multi/male/{}.paintor.log" ::: $(cut -f4 results/credible_set/PAINTOR/output_multi/male/multi_signal_loci.txt)

