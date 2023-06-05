# Rscript to remove variants that have N < median(N)

library(dplyr)
library(stringr)
library(vroom)

args = commandArgs(trailingOnly = T)

dat = vroom(args[1], col_types = cols(P = col_character()))

# If there are more than 3 studies involved in the meta-analysis, filter out
# variants that are not present in more than half of the studies, but leave
# variants that have N >= 0.2 * max(N).
#
# If there are <= 2 studies involved, keep only the variants with N >= 0.2 * max(N).
if (str_length(dat$Direction[1]) > 2) {
	min_stud = str_length(dat$Direction[1]) / 2
	dat = dat %>% mutate(N_stud = str_count(Direction, '[^\\?]')) %>% filter(N_stud >= min_stud | N >= (0.2 * max(N)))
} else {
	dat = dat %>% filter(N >= (0.2 * max(N)))
}

outname = gsub('subset', 'nfilter', args[1])

vroom_write(dat, outname)
