#! /bin/bash

mkdir -p data/neale_ukbb/neale_ldsc_sumstats/

cut -f12,16 data/neale_ukbb/ukb31063_ldsc_sumstat_manifest.tsv | grep -Fw 'TRUE' | cut -f1 | cut -d ' ' -f2 > data/neale_ukbb/neale_ukbb_url.txt

parallel -j40 "wget {} -O data/neale_ukbb/neale_ldsc_sumstats/{/.}.gz" ::: $(cat data/neale_ukbb/neale_ukbb_url.txt)

