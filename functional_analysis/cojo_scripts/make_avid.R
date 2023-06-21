# Script to make allelic value ID from summary
library(dplyr)
library(vroom)

args = commandArgs(trailingOnly = T)

dat = vroom(args[1])

minor_value = unlist(lapply(dat$minor, function(x) switch(x, 'A' = 1, 'C'= 2, 'G' = 3, 'T' = 5, -100)))
major_value = unlist(lapply(dat$major, function(x) switch(x, 'A' = 1, 'C'= 2, 'G' = 3, 'T' = 5, -100)))
value = minor_value + major_value

dat$avid = paste(dat$cpid, value, sep = '_')

dat = dat %>% select(SNP, avid)

write.table(dat, args[2], row.names = F, col.names = F, sep = '\t', quote = F)
