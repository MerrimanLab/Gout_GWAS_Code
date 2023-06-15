# Need to make sure the REF allele is actually REF allele by checking the allele with Neale's UKBB

library(dplyr)
library(vroom)

args = commandArgs(trailingOnly = T)

dat = read.table(args[1], header = T, sep = '\t', stringsAsFactors = F)

# Load in Neale's UKBB variant list with just the reference allele. Use this to
# determine if allele1 is REF or not
ref = vroom('src/slalom/ref_allele_list.tsv', col_names = c('chromosome', 'position', 'ref'))

dat = left_join(dat, ref)

# If allele1 isn't the REF allele and allele2 is the REF allele, then flip it
swap = which(dat$allele1 != dat$ref & dat$allele2 == dat$ref)

if (length(swap) > 0) {
    tmp = dat$allele2[swap]
    dat$allele2[swap] = dat$allele1[swap]
    dat$allele1[swap] = tmp
    dat$maf[swap] = 1 - dat$maf[swap]
    dat$beta[swap] = -dat$beta[swap]
    dat$z[swap] = -dat$z[swap]
}

dat = dat %>% select(rsid:p)
dat$chromosome = ifelse(dat$chromosome == 23, 'X', dat$chromosome)

# Save
filename = gsub('clean.txt', 'slalom.txt', args[1])
write.table(dat, filename, sep = '\t', col.names = T, row.names = F, quote = F)

