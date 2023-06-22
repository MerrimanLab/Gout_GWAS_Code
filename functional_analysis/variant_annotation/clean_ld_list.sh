#! /bin/bash
#################################################################################
# The PLINK output of LD is large, since it contains all tagged variants and
# not just the missense variants

parallel "cat <(head -1 results/missense/candidate_variant_proxies.{}.ld) <(grep -Fwf results/missense/missense_variants.txt results/missense/candidate_variant_proxies.{}.ld) > results/missense/candidate_variant_proxies.{}.clean.ld" ::: {eur,afr,eas,lat,tama}

