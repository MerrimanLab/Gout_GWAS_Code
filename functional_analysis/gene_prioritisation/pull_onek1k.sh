#! /bin/bash

# Pull out all eQTLs from OneK1K data

cut -f4 /Volumes/archive/merrimanlab/major_gwas_paper_archive/results_for_paper/loci_list/full_indepSNP_summary_updated_9Dec2022_withBroad.txt | tail -n+2 > data/gene_prioritisation/all_snp_list.txt
cut -f5 /Volumes/archive/merrimanlab/major_gwas_paper_archive/results_for_paper/loci_list/male_indepSNP_summary_updated_9Dec2022_withBroad.txt | tail -n+2 >> data/gene_prioritisation/all_snp_list.txt
cut -f5 /Volumes/archive/merrimanlab/major_gwas_paper_archive/results_for_paper/loci_list/female_indepSNP_summary_updated_9Dec2022_withBroad.txt | tail -n+2 >> data/gene_prioritisation/all_snp_list.txt

sort data/gene_prioritisation/all_snp_list.txt | uniq > data/gene_prioritisation/tmp && mv data/gene_prioritisation/tmp data/gene_prioritisation/all_snp_list.txt

cat <(head -1 /Volumes/archive/merrimanlab/reference_files/OneK1K/cleaned/onek1k_eqtl.memory_b_cell.tsv) <(grep -h -Fwf data/gene_prioritisation/all_snp_list.txt /Volumes/archive/merrimanlab/reference_files/OneK1K/cleaned/*) > data/gene_prioritisation/onek1k_all_snp_lead.txt

# Make a filtered version with P <= 0.05
awk -F '\t' 'NR == 1; NR > 1 && $15 <= 0.05 {print $0}' data/gene_prioritisation/onek1k_all_snp_lead.txt > data/gene_prioritisation/onek1k_all_snp_lead.p0.05.txt

# Make a filtered version with FDR <= 0.05
awk -F '\t' 'NR == 1; NR > 1 && $17 <= 0.05 {print $0}' data/gene_prioritisation/onek1k_all_snp_lead.txt > data/gene_prioritisation/onek1k_all_snp_lead.fdr0.05.txt

