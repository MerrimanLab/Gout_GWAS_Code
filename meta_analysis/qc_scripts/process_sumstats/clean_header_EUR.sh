# Very small script to rename the header of all the EUR summary data:

# CHOI
awk 'NR == 1 {print "SNP", "CHR", "BP", "minor", "major", "N", "N_case", "N_control", "MAF", "MAF_case", "MAF_control", "P", "effect", "SE"}; NR > 1 {print $0}' data/summary/EUR/choi_full_gwas_result_wgs.out | tr ' ' '\t' > data/summary/EUR/tmp && mv data/summary/EUR/tmp data/summary/EUR/choi_full_clean.out
awk 'NR == 1 {print "SNP", "CHR", "BP", "minor", "major", "N", "N_case", "N_control", "MAF", "MAF_case", "MAF_control", "P", "effect", "SE"}; NR > 1 {print $0}' data/summary/EUR/choi_female_gwas_result_wgs.out | tr ' ' '\t' > data/summary/EUR/tmp && mv data/summary/EUR/tmp data/summary/EUR/choi_female_clean.out
awk 'NR == 1 {print "SNP", "CHR", "BP", "minor", "major", "N", "N_case", "N_control", "MAF", "MAF_case", "MAF_control", "P", "effect", "SE"}; NR > 1 {print $0}' data/summary/EUR/choi_male_gwas_result_wgs.out | tr ' ' '\t' > data/summary/EUR/tmp && mv data/summary/EUR/tmp data/summary/EUR/choi_male_clean.out

# FAST/GS
awk 'NR == 1 {print "SNP", "CHR", "BP", "minor", "major", "N", "N_case", "N_control", "MAF", "MAF_case", "MAF_control", "P", "effect", "SE"}; NR > 1 {print $0}' data/summary/EUR/fast_gs_full_gwas_result_wgs.out | tr ' ' '\t' > data/summary/EUR/tmp && mv data/summary/EUR/tmp data/summary/EUR/fast_gs_full_clean.out
awk 'NR == 1 {print "SNP", "CHR", "BP", "minor", "major", "N", "N_case", "N_control", "MAF", "MAF_case", "MAF_control", "P", "effect", "SE"}; NR > 1 {print $0}' data/summary/EUR/fast_gs_male_gwas_result_wgs.out | tr ' ' '\t' > data/summary/EUR/tmp && mv data/summary/EUR/tmp data/summary/EUR/fast_gs_male_clean.out
awk 'NR == 1 {print "SNP", "CHR", "BP", "minor", "major", "N", "N_case", "N_control", "MAF", "MAF_case", "MAF_control", "P", "effect", "SE"}; NR > 1 {print $0}' data/summary/EUR/fast_gs_female_gwas_result_wgs.out | tr ' ' '\t' > data/summary/EUR/tmp && mv data/summary/EUR/tmp data/summary/EUR/fast_gs_female_clean.out

# EUROGOUT
awk 'NR == 1 {print "SNP", "CHR", "BP", "minor", "major", "N", "N_case", "N_control", "MAF", "MAF_case", "MAF_control", "P", "effect", "SE"}; NR > 1 {print $0}' data/summary/EUR/eurogout_full_gwas_result_wgs.out | tr ' ' '\t' > data/summary/EUR/tmp && mv data/summary/EUR/tmp data/summary/EUR/eurogout_full_clean.out
awk 'NR == 1 {print "SNP", "CHR", "BP", "minor", "major", "N", "N_case", "N_control", "MAF", "MAF_case", "MAF_control", "P", "effect", "SE"}; NR > 1 {print $0}' data/summary/EUR/eurogout_male_gwas_result_wgs.out | tr ' ' '\t' > data/summary/EUR/tmp && mv data/summary/EUR/tmp data/summary/EUR/eurogout_male_clean.out
awk 'NR == 1 {print "SNP", "CHR", "BP", "minor", "major", "N", "N_case", "N_control", "MAF", "MAF_case", "MAF_control", "P", "effect", "SE"}; NR > 1 {print $0}' data/summary/EUR/eurogout_female_gwas_result_wgs.out | tr ' ' '\t' > data/summary/EUR/tmp && mv data/summary/EUR/tmp data/summary/EUR/eurogout_female_clean.out

# MALMO
awk 'NR == 1 {print "SNP", "CHR", "BP", "minor", "major", "N", "N_case", "N_control", "MAF", "MAF_case", "MAF_control", "P", "effect", "SE"}; NR > 1 {print $0}' data/summary/EUR/malmo_full_gwas_result_wgs.out | tr ' ' '\t' > data/summary/EUR/tmp && mv data/summary/EUR/tmp data/summary/EUR/malmo_full_clean.out
awk 'NR == 1 {print "SNP", "CHR", "BP", "minor", "major", "N", "N_case", "N_control", "MAF", "MAF_case", "MAF_control", "P", "effect", "SE"}; NR > 1 {print $0}' data/summary/EUR/malmo_male_gwas_result_wgs.out | tr ' ' '\t' > data/summary/EUR/tmp && mv data/summary/EUR/tmp data/summary/EUR/malmo_male_clean.out
awk 'NR == 1 {print "SNP", "CHR", "BP", "minor", "major", "N", "N_case", "N_control", "MAF", "MAF_case", "MAF_control", "P", "effect", "SE"}; NR > 1 {print $0}' data/summary/EUR/malmo_female_gwas_result_wgs.out | tr ' ' '\t' > data/summary/EUR/tmp && mv data/summary/EUR/tmp data/summary/EUR/malmo_female_clean.out

# 23ANDME
awk 'NR == 1 {print "CHR", "BP", "SNP", "alleles", "major", "minor", "effect", "SE", "P", "N_control", "N_case", "MAF_control", "MAF_case", "MAF", "N" }; NR > 1 {print $0}' data/summary/EUR/europe_over17_gout_noNA_maf.dat | tr ' ' '\t' > data/summary/EUR/tmp && mv data/summary/EUR/tmp data/summary/EUR/23andMe_eur_full_clean.out
awk 'NR == 1 {print "CHR", "BP", "SNP", "alleles", "major", "minor", "effect", "SE", "P", "N_control", "N_case", "MAF_control", "MAF_case", "MAF", "N" }; NR > 1 {print $0}' data/summary/EUR/europe_female_over17_gout_noNA_maf.dat | tr ' ' '\t' > data/summary/EUR/tmp && mv data/summary/EUR/tmp data/summary/EUR/23andMe_eur_female_clean.out
awk 'NR == 1 {print "CHR", "BP", "SNP", "alleles", "major", "minor", "effect", "SE", "P", "N_control", "N_case", "MAF_control", "MAF_case", "MAF", "N" }; NR > 1 {print $0}' data/summary/EUR/europe_male_over17_gout_noNA_maf.dat | tr ' ' '\t' > data/summary/EUR/tmp && mv data/summary/EUR/tmp data/summary/EUR/23andMe_eur_male_clean.out

# KP
awk 'NR == 1 {print "CHR", "SNP", "BP", "minor", "major", "MAF", "INFO", "OR", "SE", "P", "effect", "N", "N_case", "N_control" }; NR > 1 {print $0}' data/summary/EUR/KP_full_EUR_processed.txt | tr ' ' '\t' > data/summary/EUR/tmp && mv data/summary/EUR/tmp data/summary/EUR/KP_EUR_full_clean.out
awk 'NR == 1 {print "CHR", "SNP", "BP", "minor", "major", "MAF", "INFO", "OR", "SE", "P", "effect", "N", "N_case", "N_control" }; NR > 1 {print $0}' data/summary/EUR/KP_male_EUR_processed.txt | tr ' ' '\t' > data/summary/EUR/tmp && mv data/summary/EUR/tmp data/summary/EUR/KP_EUR_male_clean.out
awk 'NR == 1 {print "CHR", "SNP", "BP", "minor", "major", "MAF", "INFO", "OR", "SE", "P", "effect", "N", "N_case", "N_control" }; NR > 1 {print $0}' data/summary/EUR/KP_female_EUR_processed.txt | tr ' ' '\t' > data/summary/EUR/tmp && mv data/summary/EUR/tmp data/summary/EUR/KP_EUR_female_clean.out

# PARTNERS
awk 'NR == 1 {print "CHR", "BP", "SNP", "minor", "major", "MAF", "effect", "SE", "P", "MAC_case", "MAC_control", "MAC", "N_case", "N_control", "N", "MAF_case", "MAF_control"}; NR > 1 {print $0}' data/summary/EUR/Partners_4sets_gout_nongout_chr1-22_MetaScore_maf.tsv | tr ' ' '\t' > data/summary/EUR/tmp && mv data/summary/EUR/tmp data/summary/EUR/Partners_full_clean.out
awk 'NR == 1 {print "CHR", "BP", "SNP", "minor", "major", "MAF", "effect", "SE", "P", "MAC_case", "MAC_control", "MAC", "N_case", "N_control", "N", "MAF_case", "MAF_control"}; NR > 1 {print $0}' data/summary/EUR/Partners_4sets_gout_nongout_M_chr1-22_MetaScore_maf.tsv | tr ' ' '\t' > data/summary/EUR/tmp && mv data/summary/EUR/tmp data/summary/EUR/Partners_male_clean.out
awk 'NR == 1 {print "CHR", "BP", "SNP", "minor", "major", "MAF", "effect", "SE", "P", "MAC_case", "MAC_control", "MAC", "N_case", "N_control", "N", "MAF_case", "MAF_control"}; NR > 1 {print $0}' data/summary/EUR/Partners_4sets_gout_nongout_F_chr1-22_MetaScore_maf.tsv | tr ' ' '\t' > data/summary/EUR/tmp && mv data/summary/EUR/tmp data/summary/EUR/Partners_female_clean.out

# UKBB
awk 'NR == 1 {print "SNP", "CHR", "BP", "minor", "major", "effect", "SE", "P", "MAF", "MAF_case", "MAF_control", "HWE", "INFO", "N", "N_case", "N_control" }; NR > 1 {print $0}' data/summary/EUR/ukbb_gwas_for_riku_study1_clean.txt | tr ' ' '\t' > data/summary/EUR/tmp && mv data/summary/EUR/tmp data/summary/EUR/ukbb_full_clean.out
awk 'NR == 1 {print "SNP", "CHR", "BP", "minor", "major", "effect", "SE", "P", "MAF", "MAF_case", "MAF_control", "HWE", "INFO", "N", "N_case", "N_control" }; NR > 1 {print $0}' data/summary/EUR/ukbb_gwas_for_riku_study2_clean.txt | tr ' ' '\t' > data/summary/EUR/tmp && mv data/summary/EUR/tmp data/summary/EUR/ukbb_male_clean.out
awk 'NR == 1 {print "SNP", "CHR", "BP", "minor", "major", "effect", "SE", "P", "MAF", "MAF_case", "MAF_control", "HWE", "INFO", "N", "N_case", "N_control" }; NR > 1 {print $0}' data/summary/EUR/ukbb_gwas_for_riku_study3_clean.txt | tr ' ' '\t' > data/summary/EUR/tmp && mv data/summary/EUR/tmp data/summary/EUR/ukbb_female_clean.out

# KOTTGEN
awk 'NR == 1 {print "CHR", "SNP", "BP", "MarkerName", "N", "minor", "major", "effect", "SE", "P", "MAF"}; NR > 1 {print $0}' data/summary/EUR/GUGC_MetaAnalysis_Results_Gout_snptracker.result_chrpos_freq.txt | tr ' ' '\t' > data/summary/EUR/tmp && mv data/summary/EUR/tmp data/summary/EUR/GUGC_full_clean.out

# FINNGEN
awk 'NR == 1 {print "CHR", "POS.b38", "major", "minor", "SNP", "gene", "P", "effect", "SE", "MAF", "MAF_case", "MAF_control", "N_case", "N_control", "N", "BP" }; NR > 1 {print $0}' data/summary/EUR/summary_stats_finngen_r3_GOUT.txt | tr ' ' '\t' > data/summary/EUR/tmp && mv data/summary/EUR/tmp data/summary/EUR/FinnGen_full_clean.out

