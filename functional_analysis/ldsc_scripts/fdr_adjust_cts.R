# Script to adjust the p-value with FDR for cell-type specific LD score
# regression results

args = commandArgs(trailingOnly = T)

dat = read.table(args[1], header = T, stringsAsFactors = F)

dat$p_fdr = p.adjust(dat$Coefficient_P_value, method = "BH")

out_name = gsub('.txt', '.fdr.txt', args[1])

write.table(dat, out_name, quote = F, sep = '\t', col.names = T, row.names = F)

