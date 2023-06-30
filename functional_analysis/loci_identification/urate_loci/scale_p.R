library(dplyr)
library(vroom)

dat = vroom('res/urate_meta/tin_ukbb.urate_meta1.clean.txt')

orig_p = vroom('res/urate_meta/original_p.txt', col_types = 'cc')
colnames(orig_p) = c('cpid', 'P_orig')

# Function to convert extreme P-values to -log10P
extreme_p = function(p) {
	if (grepl('e-', p)) {
		split_num = unlist(strsplit(p, 'e-'))
		res = -log10(as.numeric(split_num)[1]) + as.numeric(split_num[2])
	} else {
		res = -log10(as.numeric(p))
	}
	return(res)
}

orig_p$logP = unlist(lapply(orig_p$P_orig, function(x) extreme_p(x)))

# Merge this info:
dat = left_join(dat, orig_p, by = 'cpid')

# plink doesn't like extremely small P-values either, so make sure the P-values
# are just above the precision plink can understand, without changing the
# ranking of the significance of the variants (i.e. most significant variant
# should still be most significant even after adjusting the P-value):

# Convert variants with P < 1e-300 to 0:
# (Note that these P-values will/should be used ONLY for plink clumping and
# nothing else)
dat$P[which(dat$P < 1e-300)] = 0
ind = which(dat$P == 0)

# Rank these variants based on log10(P):
rank_sig = rank(-dat$logP[ind])

# The higher the rank, the smaller the P-value:
adj_p = rank_sig * .Machine$double.xmin
dat$P[ind] = adj_p

write.table(dat, 'res/urate_meta/tin_ukbb.urate_meta1.for_clumping.txt', sep = '\t', col.names = T, row.names = F, quote = F)

