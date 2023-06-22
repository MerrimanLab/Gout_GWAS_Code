#! /bin/bash

module load plink/plink1.9b6.10

# List out all the candidate causal variants (lead variants and finemapped
# variants) and pull out proxies with r2 >= 0.8

# Make a list of candidate variants from GWAS lead variants, EUR COJO analysis,
# and fine-mapping

cut -f6 /Volumes/archive/merrimanlab/major_gwas_paper_archive/results_for_paper/loci_list/full_indepSNP_summary_updated_9Dec2022_withBroad.txt | tail -n+2 > data/missense/candidate_list.txt
cut -f7 /Volumes/archive/merrimanlab/major_gwas_paper_archive/results_for_paper/loci_list/male_indepSNP_summary_updated_9Dec2022_withBroad.txt | tail -n+2 >> data/missense/candidate_list.txt
cut -f7 /Volumes/archive/merrimanlab/major_gwas_paper_archive/results_for_paper/loci_list/female_indepSNP_summary_updated_9Dec2022_withBroad.txt | tail -n+2 >> data/missense/candidate_list.txt
cut -f3 results/credible_set/combined_results/top_candidate_cred_var.combined.txt | tail -n+2 >> data/missense/candidate_list.txt

# Get COJO variants
Rscript src/missense/get_cojo_var.R data/cojo_dat/full_loci_summary_1jun2022_cojo.txt

# Combine and clean the list
cat data/missense/candidate_list.txt data/missense/cojo_vars.txt | sort | uniq > data/missense/tmp && mv data/missense/tmp data/missense/candidate_list.txt

# Generate proxies using 1KGP data
plink --bfile data/1kgp_ref/EUR_wgs --show-tags data/missense/candidate_list.txt --list-all --tag-r2 0.8 --tag-kb 500 --out results/missense/candidate_variant_proxies.eur
plink --bfile data/1kgp_ref/AFR_wgs --show-tags data/missense/candidate_list.txt --list-all --tag-r2 0.8 --tag-kb 500 --out results/missense/candidate_variant_proxies.afr
plink --bfile data/1kgp_ref/EAS_wgs --show-tags data/missense/candidate_list.txt --list-all --tag-r2 0.8 --tag-kb 500 --out results/missense/candidate_variant_proxies.eas
plink --bfile data/1kgp_ref/LAT_wgs --show-tags data/missense/candidate_list.txt --list-all --tag-r2 0.8 --tag-kb 500 --out results/missense/candidate_variant_proxies.lat
plink --bfile data/1kgp_ref/TAMA_wgs --show-tags data/missense/candidate_list.txt --list-all --tag-r2 0.8 --tag-kb 500 --out results/missense/candidate_variant_proxies.tama

# Make a list of all the tag variants:
cat results/missense/candidate_variant_proxies.*.tags | sort | uniq > results/missense/candidate_variant_proxies.all.txt

plink --bfile data/1kgp_ref/EUR_wgs  --r2 --ld-window 50000 --ld-window-kb 500 --ld-window-r2 0 --ld-snp-list results/missense/candidate_variant_proxies.all.txt --out results/missense/candidate_variant_proxies.eur
plink --bfile data/1kgp_ref/AFR_wgs  --r2 --ld-window 50000 --ld-window-kb 500 --ld-window-r2 0 --ld-snp-list results/missense/candidate_variant_proxies.all.txt --out results/missense/candidate_variant_proxies.afr
plink --bfile data/1kgp_ref/EAS_wgs  --r2 --ld-window 50000 --ld-window-kb 500 --ld-window-r2 0 --ld-snp-list results/missense/candidate_variant_proxies.all.txt --out results/missense/candidate_variant_proxies.eas
plink --bfile data/1kgp_ref/LAT_wgs  --r2 --ld-window 50000 --ld-window-kb 500 --ld-window-r2 0 --ld-snp-list results/missense/candidate_variant_proxies.all.txt --out results/missense/candidate_variant_proxies.lat
plink --bfile data/1kgp_ref/TAMA_wgs --r2 --ld-window 50000 --ld-window-kb 500 --ld-window-r2 0 --ld-snp-list results/missense/candidate_variant_proxies.all.txt --out results/missense/candidate_variant_proxies.tama
