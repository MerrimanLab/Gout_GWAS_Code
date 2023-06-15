# Script to generate correlation matrix of the significant annotations from
# sLDSC and output the annotations to be used for PAINTOR

library(dplyr)
library(purrr)

sig_annot = read.table('results/ldsc/paintor_ldsc/sig_annot.txt', sep = '\t', header = T, stringsAsFactors = F)

files = paste('data/ldsc/1kgp_ref_paintor/tmp_annot', sig_annot$annot_origin, sig_annot$annotation, sep = '/')

annot_list = list()

for (i in 1:length(files)) {
	prefix = files[i]
	file_name = paste(prefix, 'chr', 1:22, '.annot.gz', sep = '')
	annot = map_dfr(file_name, ~ read.table(gzfile(.x), header = T, sep = '\t', stringsAsFactors = F))
	annot_list[[i]] = annot
}

annot_mat = map_dfc(annot_list, ~ .x)
colnames(annot_mat) = files

annot_cor = cor(annot_mat, method = 'pearson') ^ 2
diag(annot_cor) = 0

write.table(annot_cor, 'results/ldsc/paintor_ldsc/annotation_correlation.txt', col.names = T, row.names = F, quote = F)

# Pull out a set of annotations that are uncorrelated to one another (squared
# correlation <= 0.2)
count = 1
for (i in 1:nrow(annot_cor) - 1) {
    tmp = annot_cor[1:count, 1:count]
    if (any(tmp >= 0.2)) {
        annot_cor = annot_cor[-count, -count]
        count = count - 1
    }
    count = count + 1
}

# Save list of uncorrelated annotations:
writeLines(colnames(annot_cor), 'results/ldsc/paintor_ldsc/sig_annot.uncorrelated.txt')
