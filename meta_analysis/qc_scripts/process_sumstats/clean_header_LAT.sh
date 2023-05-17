# Very small script to rename the header of all the LAT summary data:

# 23ANDME
awk 'NR == 1 {print "CHR", "BP", "SNP", "alleles", "major", "minor", "effect", "SE", "P", "N_control", "N_case", "MAF_control", "MAF_case", "MAF", "N" }; NR > 1 {print $0}' data/summary/LAT/latino_over17_gout_noNA_maf.dat | tr ' ' '\t' > data/summary/LAT/tmp && mv data/summary/LAT/tmp data/summary/LAT/23andMe_lat_full_clean.out
awk 'NR == 1 {print "CHR", "BP", "SNP", "alleles", "major", "minor", "effect", "SE", "P", "N_control", "N_case", "MAF_control", "MAF_case", "MAF", "N" }; NR > 1 {print $0}' data/summary/LAT/latino_male_over17_gout_noNA_maf.dat | tr ' ' '\t' > data/summary/LAT/tmp && mv data/summary/LAT/tmp data/summary/LAT/23andMe_lat_male_clean.out
awk 'NR == 1 {print "CHR", "BP", "SNP", "alleles", "major", "minor", "effect", "SE", "P", "N_control", "N_case", "MAF_control", "MAF_case", "MAF", "N" }; NR > 1 {print $0}' data/summary/LAT/latino_female_over17_gout_noNA_maf.dat | tr ' ' '\t' > data/summary/LAT/tmp && mv data/summary/LAT/tmp data/summary/LAT/23andMe_lat_female_clean.out

# KP
awk 'NR == 1 {print "CHR", "SNP", "BP", "minor", "major", "MAF", "INFO", "OR", "SE", "P", "effect", "N", "N_case", "N_control"}; NR > 1 {print $0}' data/summary/LAT/KP_full_LAT_processed.txt | tr ' ' '\t' > data/summary/LAT/tmp && mv data/summary/LAT/tmp data/summary/LAT/KP_LAT_full_clean.out
awk 'NR == 1 {print "CHR", "SNP", "BP", "minor", "major", "MAF", "INFO", "OR", "SE", "P", "effect", "N", "N_case", "N_control"}; NR > 1 {print $0}' data/summary/LAT/KP_female_LAT_processed.txt | tr ' ' '\t' > data/summary/LAT/tmp && mv data/summary/LAT/tmp data/summary/LAT/KP_LAT_female_clean.out
awk 'NR == 1 {print "CHR", "SNP", "BP", "minor", "major", "MAF", "INFO", "OR", "SE", "P", "effect", "N", "N_case", "N_control"}; NR > 1 {print $0}' data/summary/LAT/KP_male_LAT_processed.txt | tr ' ' '\t' > data/summary/LAT/tmp && mv data/summary/LAT/tmp data/summary/LAT/KP_LAT_male_clean.out
