# Short script to generate manhattan and QQ plots using METAL output.
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
data$P = as.numeric(data[,'P-value'])
data = data[!is.na(data$P), ]

out_name = gsub('\\..*', '', args)

png(paste(out_name, 'QQ.png', sep = '_'), width = 1000, height = 1000)
qq(data$P, main = "QQ plot")
dev.off()

png(paste(out_name, 'manhattan.png', sep = '_'), width = 2000, height = 1000)
manhattan(data, chr = "CHR", bp = "BP", p = "P", snp = "MarkerName", main = "Manhattan plot")
dev.off()

sub_data = data[which(data[,"P"] >= 1e-30), ]

png(paste(out_name, 'manhattan_zoom.png', sep = '_'), width = 2000, height = 1000)
manhattan(sub_data, chr = "CHR", bp = "BP", p = "P", snp = "MarkerName", main = "Manhattan plot")
dev.off()
