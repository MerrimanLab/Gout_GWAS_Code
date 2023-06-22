# Script to merge all the PPH4 from all the coloc results for gout vs urate
# colocalisation
library(dplyr)
library(tidyr)
library(vroom)
library(purrr)

# Load coloc results
files = list.files('res/coloc/gout_urate/', pattern = 'coloc_res.txt', full.names = T)
res = map_dfr(files, ~ read.table(.x, sep = '\t', header = T, stringsAsFactors = F))

# Load loci info
loci = read.table('res/urate_clumping/tin_ukbb.loci_list.txt', sep = '\t', header = T, stringsAsFactors = F)

# Add CHR/POS info
urate = vroom('res/urate_meta/tin_ukbb.urate_meta1.clean.txt')

loci = urate %>% select(SNP, CHR, POS) %>% left_join(loci, ., by = c('lead' = 'SNP'))
loci = loci %>% select(loci, lead, CHR, POS) %>% arrange(CHR, POS)
colnames(loci)[2] = 'SNP'

# Join the loci info and remove any that have NA values in CHR/POS, since these
# will be from the gout loci and not the urate loci
res = right_join(loci, res) %>% filter(!is.na(CHR))

# Filter out for loci only in urate (i.e. PP.H2.abf >= 0.8)
res_urate = res %>% filter(PP.H2.abf >= 0.8)

# What about those that may have different signals?
res_diff = res %>% filter(PP.H3.abf >= 0.8)

# Also pull out those that have H2 + H4 >= 0.8 - there may be some loci that
# have significant signal at urate and sub-significant signal for gout, meaning
# that the locus may have distributed the probability to H2 and H4, where
# neither passes the 0.8 threshold
res_weak = res %>% filter(PP.H2.abf < 0.8 & PP.H4.abf < 0.8) %>% mutate(h2_plus_h4 = PP.H2.abf + PP.H4.abf) %>% filter(h2_plus_h4 >= 0.8)

# Save results
write.table(res, 'res/coloc/gout_urate/coloc_res.urate_loci.txt', sep = '\t', col.names = T, row.names = F, quote = F)
write.table(res_urate, 'res/coloc/gout_urate/coloc_res.h2.urate_loci.txt', sep = '\t', col.names = T, row.names = F, quote = F)
write.table(res_diff, 'res/coloc/gout_urate/coloc_res.h3.urate_loci.txt', sep = '\t', col.names = T, row.names = F, quote = F)
write.table(res_weak, 'res/coloc/gout_urate/coloc_res.h2_plus_h4.urate_loci.txt', sep = '\t', col.names = T, row.names = F, quote = F)
