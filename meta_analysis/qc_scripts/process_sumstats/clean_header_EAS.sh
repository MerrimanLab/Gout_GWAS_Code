# Very small script to rename the header of all the EAS summary data:

# KOREAN:
awk 'NR == 1 {print "CHR", "SNP", "BP", "minor", "major", "MAF", "MAF_case", "MAF_control", "CHISQ", "P", "OR", "SE", "effect", "N", "N_case", "N_control" }; NR > 1 {print $0}' data/summary/EAS/korean_gout.txt | tr ' ' '\t' > data/summary/EAS/tmp && mv data/summary/EAS/tmp data/summary/EAS/korean_full_clean.out

# 23ANDME:
awk 'NR == 1 {print "CHR", "BP", "SNP", "alleles", "major", "minor", "effect", "SE", "P", "N_control", "N_case", "MAF_control", "MAF_case", "MAF", "N" }; NR > 1 {print $0}' data/summary/EAS/east_asian_over17_gout_noNA_maf.dat | tr ' ' '\t' > data/summary/EAS/tmp && mv data/summary/EAS/tmp data/summary/EAS/23andMe_eas_full_clean.out
awk 'NR == 1 {print "CHR", "BP", "SNP", "alleles", "major", "minor", "effect", "SE", "P", "N_control", "N_case", "MAF_control", "MAF_case", "MAF", "N" }; NR > 1 {print $0}' data/summary/EAS/east_asian_male_over17_gout_noNA_maf.dat | tr ' ' '\t' > data/summary/EAS/tmp && mv data/summary/EAS/tmp data/summary/EAS/23andMe_eas_male_clean.out
awk 'NR == 1 {print "CHR", "BP", "SNP", "alleles", "major", "minor", "effect", "SE", "P", "N_control", "N_case", "MAF_control", "MAF_case", "MAF", "N" }; NR > 1 {print $0}' data/summary/EAS/east_asian_female_over17_gout_noNA_maf.dat | tr ' ' '\t' > data/summary/EAS/tmp && mv data/summary/EAS/tmp data/summary/EAS/23andMe_eas_female_clean.out

# JAPANESE
awk 'NR == 1 {print "SNP", "CHR", "BP", "minor", "major", "Rsq", "MAF_case", "MAF_control", "effect", "SE", "P", "N_case", "N_control", "N", "MAC_case", "MAC_control", "MAC_total", "MAF"}; NR > 1 {print $0}' data/summary/EAS/jap_japonica_maf.tsv | tr ' ' '\t' > data/summary/EAS/tmp && mv data/summary/EAS/tmp data/summary/EAS/japonica_full_clean.out
awk 'NR == 1 {print "SNP", "CHR", "BP", "minor", "major", "Rsq", "MAF_case", "MAF_control", "effect", "SE", "P", "N_case", "N_control", "N", "MAC_case", "MAC_control", "MAC_total", "MAF"}; NR > 1 {print $0}' data/summary/EAS/jap_illumina_maf.tsv | tr ' ' '\t' > data/summary/EAS/tmp && mv data/summary/EAS/tmp data/summary/EAS/illumina_full_clean.out
cp data/summary/EAS/japonica_full_clean.out data/summary/EAS/japonica_male_clean.out
cp data/summary/EAS/illumina_full_clean.out data/summary/EAS/illumina_male_clean.out

# CHINESE
awk 'NR == 1 {print "CHR", "SNP", "BP", "major", "minor", "P", "effect", "SE", "MAF_case", "MAF_control", "MAF", "N"}; NR > 1 {print $0}' data/summary/EAS/chinese_gout.tsv | tr ' ' '\t' > data/summary/EAS/tmp && mv data/summary/EAS/tmp data/summary/EAS/chinese_full_clean.out
awk 'NR == 1 {print "CHR", "SNP", "BP", "minor", "major", "MAF_case", "MAF_control", "OR", "SE", "P", "A1_A2", "MAC_case", "MAC_control", "MAC_total", "N", "N_case", "N_control", "MAF", "effect"}; NR > 1 {print $0}' data/summary/EAS/chinese_gout_male.tsv | tr ' ' '\t' > data/summary/EAS/tmp && mv data/summary/EAS/tmp data/summary/EAS/chinese_male_clean.out

# KP
awk 'NR == 1 {print "CHR", "SNP", "BP", "minor", "major", "MAF", "INFO", "OR", "SE", "P", "effect", "N", "N_case", "N_control" }; NR > 1 {print $0}' data/summary/EAS/KP_full_EAS_processed.txt | tr ' ' '\t' > data/summary/EAS/tmp && mv data/summary/EAS/tmp data/summary/EAS/KP_EAS_full_clean.out
awk 'NR == 1 {print "CHR", "SNP", "BP", "minor", "major", "MAF", "INFO", "OR", "SE", "P", "effect", "N", "N_case", "N_control" }; NR > 1 {print $0}' data/summary/EAS/KP_male_EAS_processed.txt | tr ' ' '\t' > data/summary/EAS/tmp && mv data/summary/EAS/tmp data/summary/EAS/KP_EAS_male_clean.out
awk 'NR == 1 {print "CHR", "SNP", "BP", "minor", "major", "MAF", "INFO", "OR", "SE", "P", "effect", "N", "N_case", "N_control" }; NR > 1 {print $0}' data/summary/EAS/KP_female_EAS_processed.txt | tr ' ' '\t' > data/summary/EAS/tmp && mv data/summary/EAS/tmp data/summary/EAS/KP_EAS_female_clean.out
