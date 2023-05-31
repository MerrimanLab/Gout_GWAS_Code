# Short script to generate manhattan and QQ plots using log10P-converted METAL
# output.
#
# Note that the SNP/variant ID should be in CHR_POS format.

library(vroom)
library(dplyr)
library(qqman)

args = commandArgs(trailingOnly = T)

data = as.data.frame(vroom(args[1]))

chr_bp_inf = data$MarkerName %>% strsplit('_') %>% unlist %>% as.numeric %>% matrix(., ncol = 2, byrow = T)

data$CHR = chr_bp_inf[,1]
data$BP = chr_bp_inf[,2]
data = data[!is.na(data$log10P), ]

out_name = gsub('\\..*', '', args)

png(paste(out_name, 'QQ.png', sep = '_'), width = 1000, height = 1000)
qq(as.numeric(data[,'P-value']), main = "QQ plot")
dev.off()

sub_data = data[which(data[,"log10P"] <= 30), ]

png(paste(out_name, 'manhattan_zoom.png', sep = '_'), width = 2000, height = 1000)
manhattan(sub_data, chr = "CHR", bp = "BP", p = "log10P", snp = "MarkerName", main = "Manhattan plot", logp = F)
dev.off()

png(paste(out_name, 'manhattan.png', sep = '_'), width = 2000, height = 1000)
manhattan(data, chr = "CHR", bp = "BP", p = "log10P", snp = "MarkerName", main = "Manhattan plot", logp = F)
dev.off()

