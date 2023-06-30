#! /bin/bash
# Run LD of all the clump lead variants

module load plink/plink1.9b6.10

cut -f3 res/urate_clumping/all_clump_lead.txt | tail -n+2 | sort | uniq > res/urate_clumping/clump_lead_ld.input.txt

plink --bfile dat/1kgp/EUR_wgs.no_relatives.no_indel.biallelic --extract res/urate_clumping/clump_lead_ld.input.txt --ld-snp-list res/urate_clumping/clump_lead_ld.input.txt --r2 inter-chr --ld-window-r2 0 --out res/urate_clumping/clump_lead_ld

