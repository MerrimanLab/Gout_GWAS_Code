# Process the SLALOM result

library(dplyr)
library(purrr)
library(vroom)

args = commandArgs(trailingOnly = T)

file_list = list.files(args[1], pattern = 'result.txt', full.names = T)
dat_list = map(file_list, ~ read.table(.x, sep = '\t', header = T, stringsAsFactors = F))

# Name the list with the corresponding locus name/rsid:
locus = gsub('.*/', '', file_list)
locus = gsub('\\..*', '', locus)
names(dat_list) = locus

# Function to process the SLALOM data and determine if a locus is suspicious or
# not
check_susp = function(dat, r2 = 0.6, log_dent = 4) {
    # Remove variants with no correlation with the lead variant (i.e. not
    # present in GnomAD)
    dat = dat[!is.na(dat$r),]
    # Remove variants that are in 1.0 r2 with the lead variant, since this
    # ruins the DENTIST-S score (makes it Inf due to the denominator being 0)
    dat = dat[dat$r != 1,]
    dat$r2 = dat$r ^ 2
    # Check if any variant has r2 > 0.6 and -log10(DENTIST-S p-value) > 4:
    susp = any(dat$r2 > r2 & dat$nlog10p_dentist_s > log_dent)
    return(susp)
}

susp = map_lgl(dat_list, ~ check_susp(.x))

# For checking how many loci were suspicious:
# which(susp) %>% length
# length(dat_list)

# Make a function to pull out variants that had 100% LD with the lead variant
# and/or violated the SLALOM threshold (r2 > 0.6 and -log10(DENTIST p) > 4)
get_susp_var = function(dat, r2 = 0.6, log_dent = 4) {
    dat = dat[!is.na(dat$r),]
    dat$r2 = dat$r ^ 2
    ind = unique(c(which(dat$r == 1), which(dat$r2 > r2 & dat$nlog10p_dentist_s > log_dent)))
    if(length(ind) > 0) {
        res = dat[ind,]
    }
    return(res)
}

# Pull out variants with 100% LD with lead and violating variants from the
# suspicious loci
susp_loci = dat_list[which(susp)]
susp_var = map(susp_loci, ~ get_susp_var(.x))

# Save the list of suspicious variants in these loci
for (i in 1:length(susp_var)) {
    filename = paste(c(args[1], names(susp_var)[i], ".suspicious_variants.txt"), collapse = "")
    write.table(susp_var[[i]], filename, sep = '\t', col.names = T, row.names = F, quote = F)
}

# Also save the list of loci that were suspicious
writeLines(names(susp_loci), paste(args[1], 'suspicious_loci.txt', sep = '/'))
