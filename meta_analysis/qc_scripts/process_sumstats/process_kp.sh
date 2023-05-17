#! /bin/bash
# Small script to organise the KP data sets.

# Add sample size, calculate beta, and change chr25 to chr23 for KP data:

# AFR
awk 'NR == 1 {print $0, "effect", "N", "N_case", "N_control"}; NR > 1 {print $0, log($8), 2728, 258, 2470}' data/summary/AFR/Gout_AllControls_AFRImputed_allnoNA.txt | sed 's/^25/23/g' > data/summary/AFR/KP_full_AFR_processed.txt
awk 'NR == 1 {print $0, "effect", "N", "N_case", "N_control"}; NR > 1 {print $0, log($8), 1573, 88,1485}' data/summary/AFR/WomenGout_AllControls_AFRImputed_allnoNA.txt | sed 's/^25/23/g' > data/summary/AFR/KP_female_AFR_processed.txt
awk 'NR == 1 {print $0, "effect", "N", "N_case", "N_control"}; NR > 1 {print $0, log($8), 1155, 170, 985}' data/summary/AFR/MenGout_AllControls_AFRImputed_allnoNA.txt | sed 's/^25/23/g' > data/summary/AFR/KP_male_AFR_processed.txt

# EAS
awk 'NR == 1 {print $0, "effect", "N", "N_case", "N_control"}; NR > 1 {print $0, log($8), 6067, 647, 5420}' data/summary/EAS/Gout_AllControls_EASImputed_allnoNA.txt | sed 's/^25/23/g' > data/summary/EAS/KP_full_EAS_processed.txt
awk 'NR == 1 {print $0, "effect", "N", "N_case", "N_control"}; NR > 1 {print $0, log($8), 3410, 125, 3285}' data/summary/EAS/WomenGout_AllControls_EASImputed_allnoNA.txt | sed 's/^25/23/g' > data/summary/EAS/KP_female_EAS_processed.txt
awk 'NR == 1 {print $0, "effect", "N", "N_case", "N_control"}; NR > 1 {print $0, log($8), 2657, 522, 2135}' data/summary/EAS/MenGout_AllControls_EASImputed_allnoNA.txt | sed 's/^25/23/g' > data/summary/EAS/KP_male_EAS_processed.txt

# EUR
awk 'NR == 1 {print $0, "effect", "N", "N_case", "N_control"}; NR > 1 {print $0, log($8), 69042, 4232, 64810}' data/summary/EUR/Gout_AllControls_EURImputed_allnoNA.txt | sed 's/^25/23/g' > data/summary/EUR/KP_full_EUR_processed.txt
awk 'NR == 1 {print $0, "effect", "N", "N_case", "N_control"}; NR > 1 {print $0, log($8), 29799, 3162, 26637}' data/summary/EUR/MenGout_AllControls_EURImputed_allnoNA.txt | sed 's/^25/23/g' > data/summary/EUR/KP_male_EUR_processed.txt
awk 'NR == 1 {print $0, "effect", "N", "N_case", "N_control"}; NR > 1 {print $0, log($8), 39243, 1070, 38173}' data/summary/EUR/WomenGout_AllControls_EURImputed_allnoNA.txt | sed 's/^25/23/g' > data/summary/EUR/KP_female_EUR_processed.txt

# LAT
awk 'NR == 1 {print $0, "effect", "N", "N_case", "N_control"}; NR > 1 {print $0, log($8), 5742, 268, 5474}' data/summary/LAT/Gout_AllControls_LATImputed_allnoNA.txt | sed 's/^25/23/g' > data/summary/LAT/KP_full_LAT_processed.txt
awk 'NR == 1 {print $0, "effect", "N", "N_case", "N_control"}; NR > 1 {print $0, log($8), 3425, 68, 3357}' data/summary/LAT/WomenGout_AllControls_LATImputed_allnoNA.txt | sed 's/^25/23/g' > data/summary/LAT/KP_female_LAT_processed.txt
awk 'NR == 1 {print $0, "effect", "N", "N_case", "N_control"}; NR > 1 {print $0, log($8), 2317, 200, 2117}' data/summary/LAT/MenGout_AllControls_LATImputed_allnoNA.txt | sed 's/^25/23/g' > data/summary/LAT/KP_male_LAT_processed.txt
