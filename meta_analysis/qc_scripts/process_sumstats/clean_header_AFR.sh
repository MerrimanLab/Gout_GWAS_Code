# Very small script to rename the header of all the AFR summary data:

# 23ANDME
awk 'NR == 1 {print "CHR", "BP", "SNP", "alleles", "major", "minor", "effect", "SE", "P", "N_control", "N_case", "MAF_control", "MAF_case", "MAF", "N" }; NR > 1 {print $0}' data/summary/AFR/african_american_over17_gout_noNA_maf.dat | tr ' ' '\t' > data/summary/AFR/tmp && mv data/summary/AFR/tmp data/summary/AFR/23andMe_afr_full_clean.out
awk 'NR == 1 {print "CHR", "BP", "SNP", "alleles", "major", "minor", "effect", "SE", "P", "N_control", "N_case", "MAF_control", "MAF_case", "MAF", "N" }; NR > 1 {print $0}' data/summary/AFR/african_american_male_over17_gout_noNA_maf.dat | tr ' ' '\t' > data/summary/AFR/tmp && mv data/summary/AFR/tmp data/summary/AFR/23andMe_afr_male_clean.out
awk 'NR == 1 {print "CHR", "BP", "SNP", "alleles", "major", "minor", "effect", "SE", "P", "N_control", "N_case", "MAF_control", "MAF_case", "MAF", "N" }; NR > 1 {print $0}' data/summary/AFR/african_american_female_over17_gout_noNA_maf.dat | tr ' ' '\t' > data/summary/AFR/tmp && mv data/summary/AFR/tmp data/summary/AFR/23andMe_afr_female_clean.out

# KP
awk 'NR == 1 {print "CHR", "SNP", "BP", "minor", "major", "MAF", "INFO", "OR", "SE", "P", "effect", "N", "N_case", "N_control"}; NR > 1 {print $0}' data/summary/AFR/KP_full_AFR_processed.txt | tr ' ' '\t' > data/summary/AFR/tmp && mv data/summary/AFR/tmp data/summary/AFR/KP_AFR_full_clean.out
awk 'NR == 1 {print "CHR", "SNP", "BP", "minor", "major", "MAF", "INFO", "OR", "SE", "P", "effect", "N", "N_case", "N_control"}; NR > 1 {print $0}' data/summary/AFR/KP_female_AFR_processed.txt | tr ' ' '\t' > data/summary/AFR/tmp && mv data/summary/AFR/tmp data/summary/AFR/KP_AFR_female_clean.out
awk 'NR == 1 {print "CHR", "SNP", "BP", "minor", "major", "MAF", "INFO", "OR", "SE", "P", "effect", "N", "N_case", "N_control"}; NR > 1 {print $0}' data/summary/AFR/KP_male_AFR_processed.txt | tr ' ' '\t' > data/summary/AFR/tmp && mv data/summary/AFR/tmp data/summary/AFR/KP_AFR_male_clean.out
