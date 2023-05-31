################################################################################
# Inter-study comparison stuff
################################################################################

library(vroom)
library(dplyr)
library(tidyr)
library(purrr)
#devtools::install_github('DavisVaughan/furrr')
library(furrr)
library(qqman)

# Setup multicore:
plan(strategy = 'multicore', workers = 10)
options(future.globals.maxSize= 50 * 1024 ^ 3 )

args = commandArgs(trailingOnly = T)

data_list = list()

dat_names = gsub('.*/', '', args)
dat_names = gsub('(^[^_]*)(.*$)', '\\1', dat_names)

# Load data:
for (i in 1:length(args)) {
	data_list[[i]] = vroom(args[i])
}

names(data_list) = dat_names

reorder = order(unlist(lapply(data_list, nrow)), decreasing = T)

data_list = lapply(data_list[reorder], as.data.frame)
dat_names = dat_names[reorder]

# Function to compare the parameters of two studies passed to the function:
compare_studies = function(dat1, dat2, names) {
	# Find common variants between two studies
	id1 = dat1$cpid[which(dat1$cpid %in% dat2$cpid)]
	id2 = dat2$cpid[which(dat2$cpid %in% dat1$cpid)]
	common_snps = id1[which(id1 %in% id2)]

	# Subset and order everything based on first data:
	dat1 = dat1[dat1$cpid %in% common_snps, ]
	dat2 = dat2[dat2$cpid %in% common_snps, ]
	rownames(dat1) = dat1$cpid
	rownames(dat2) = dat2$cpid
	dat2 = dat2[rownames(dat1),]

	# Generate A1_A2 column:
	dat1$A1_A2 = paste(dat1$minor, dat1$major, sep = '_')
	dat2$A1_A2 = paste(dat2$minor, dat2$major, sep = '_')

	# Switch MAF, effects, and alleles based on dat1, so they are consistent
	# between the studies:
	ind = which(dat1$A1_A2 != dat2$A1_A2)

	tmp = dat2[ind, "minor"]
	dat2[ind, "minor"] = dat2[ind, "major"]
	dat2[ind, "major"] = tmp
	dat2[ind, "effect"] = -dat2[ind, "effect"]
	dat2[ind, "MAF"] = 1 - dat2[ind, "MAF"]
	dat2$A1_A2 = paste(dat2[,"minor"], dat2[,"major"], sep = '_')

	# Setup filenames:
	tmp_names = gsub('.*/', '', names)
	plot_basename = paste(tmp_names[1], tmp_names[2], sep = " vs. ")
	file_basename = paste(names[1], tmp_names[2], sep = "_vs_")

	# Deal with case/control MAF plots only if both data sets contain
	# case/control MAF info:
	if (any(grepl('MAF_case|MAF_control', colnames(dat1))) & any(grepl('MAF_case|MAF_control', colnames(dat2)))) {
		dat2[ind, "MAF_case"] = 1 - dat2[ind, "MAF_case"]
		dat2[ind, "MAF_control"] = 1 - dat2[ind, "MAF_control"]

		# temporary data for plotting
		tmp_case = bind_cols(MAF1 = dat1$MAF_case, MAF2 = dat2$MAF_case) %>% drop_na() %>% distinct
		tmp_control = bind_cols(MAF1 = dat1$MAF_control, MAF2 = dat2$MAF_control) %>% drop_na() %>% distinct

		# Case/control MAF plots:
		png(paste(file_basename, 'MAF_case.png', sep = '_'), width = 1000, height = 1000)
		plot(x = tmp_case$MAF1, y = tmp_case$MAF2, main = paste(plot_basename, "MAF (case)", sep = ' '), xlab = names[1], ylab = names[2], pch = 20)
		abline(a = 0, b = 1, col = "red")
		dev.off()

		png(paste(file_basename, 'MAF_control.png', sep = '_'), width = 1000, height = 1000)
		plot(x = tmp_control$MAF1, y = tmp_control$MAF2, main = paste(plot_basename, "MAF (control)", sep = ' '), xlab = names[1], ylab = names[2], pch = 20)
		abline(a = 0, b = 1, col = "red")
		dev.off()
	}

	# Generate scatter plots for the important parameters of two studies
	png(paste(file_basename, 'logOR.png', sep = '_'), width = 1000, height = 1000)
	tmp_effect = bind_cols(effect1 = dat1$effect, effect2 = dat2$effect) %>% drop_na() %>% distinct
	plot(x = tmp_effect$effect1, y = tmp_effect$effect2, main = paste(plot_basename, "logOR", sep = ' '), xlab = names[1], ylab = names[2], pch = 20)
	abline(h = 0, v = 0, col = "red")
	dev.off()

	lowf_ind = which(!(dat1$MAF < 0.01 | dat1$MAF > 0.99 | dat2$MAF < 0.01 | dat2$MAF > 0.99))

	png(paste(file_basename, 'logOR_noLF.png', sep = '_'), width = 1000, height = 1000)
	tmp_effect = bind_cols(effect1 = dat1$effect[lowf_ind], effect2 = dat2$effect[lowf_ind]) %>% drop_na() %>% distinct
	plot(x = tmp_effect$effect1, y = tmp_effect$effect2, main = paste(plot_basename, "logOR", sep = ' '), xlab = names[1], ylab = names[2], pch = 20)
	abline(h = 0, v = 0, col = "red")
	dev.off()

	png(paste(file_basename, 'MAF.png', sep = '_'), width = 1000, height = 1000)
	tmp = bind_cols(MAF1 = dat1$MAF, MAF2 = dat2$MAF) %>% drop_na() %>% distinct
	plot(x = tmp$MAF1, y = tmp$MAF2, main = paste(plot_basename, "MAF", sep = ' '), xlab = names[1], ylab = names[2], pch = 20)
	abline(a = 0, b = 1, col = "red")
	dev.off()

}

# Compare studies:
combination = t(combn(c(1:length(data_list)), 2))

# Generate vector of output name:
ancestry = gsub('\\/[^\\/]*$', '', args)
ancestry = unique(gsub('.*/', '', ancestry))

sample_group = ifelse(grepl("full", args[1]), "full/", ifelse(grepl("female", args[1]), "female/", "male/"))

out_dir = paste('results/pre_meta_qc', ancestry, sep = '/')
out_dir = paste(out_dir, sample_group, sep = '/')

out_name = paste(out_dir, dat_names, sep = '')

future_walk2(combination[, 1], combination[, 2], ~ compare_studies(data_list[[.x]], data_list[[.y]], out_name[c(.x, .y)]))

# Function to generate distribution plots for a given study:
plot_dist = function(dat, data_name, out_name) {
	# Plot beta/logOR
	png(paste(out_name, 'logOR_distribution.png', sep = '_'), width = 1000, height = 1000)
	plot(density(dat$effect), main = paste("Density plot of logOR for", data_name), xlab = "logOR")
	dev.off()

	# Plot P-Z
	png(paste(out_name, 'PvsZ.png', sep = '_'), width = 1000, height = 1000)
	plot(x = -log10(dat$P), y = -log10(2 * pnorm(-abs(dat$effect/dat$SE))), main = paste("P vs. Z plot for", data_name), xlab = "-log10(P)", ylab = "-log10(P.Ztest)")
	abline(0, 1, col = "red")
	dev.off()

	# Plot some data without low frequency variants (MAF > 0.01)
	tmp_dat = dat[dat$MAF < 0.99 & dat$MAF > 0.01, ]

	png(paste(out_name, 'PvsZ_noLF.png', sep = '_'), width = 1000, height = 1000)
	plot(x = -log10(tmp_dat$P), y = -log10(2 * pnorm(-abs(tmp_dat$effect/tmp_dat$SE))), main = paste("P vs. Z plot for", data_name), xlab = "-log10(P)", ylab = "-log10(P.Ztest)")
	abline(0, 1, col = "red")
	dev.off()

	# Generate QQ and manhattan plots:
	png(paste(out_name, 'QQ.png', sep = '_'), width = 1000, height = 1000)
	qq(dat$P, main = paste("QQ plot of", data_name))
	dev.off()

	# Reduce number of SNPs for manhattan plot:
	small_dat = dat[dat$P < 0.2, ]
	small_dat$P[small_dat$P == 0] = .Machine$double.xmin

	png(paste(out_name, 'manhattan.png', sep = '_'), width = 2000, height = 1000)
	small_dat$CHR = as.integer(small_dat$CHR)
	small_dat$BP = as.integer(small_dat$BP)
	manhattan(small_dat, chr = "CHR", bp = "BP", p = "P", snp = "SNP", main = paste("Manhattan plot of", data_name))
	dev.off()
}

out_name = gsub('_clean.*', '', args[reorder])
out_name = gsub('.*/', '', out_name)
out_name = paste(out_dir, out_name, sep = '')
names(out_name) = names(data_list)

future_iwalk(data_list, ~ plot_dist(.x, .y, as.vector(out_name[.y])))

# Compare quality of summary stats:

# First make a list of cpid common across all data sets:
common_snps = data_list[[1]]$cpid

for (i in 2:length(data_list)) {
	id1 = data_list[[i]]$cpid[which(data_list[[i]]$cpid %in% common_snps)]
	id2 = common_snps[which(common_snps %in% data_list[[i]]$cpid)]
	common_snps = id1[which(id1 %in% id2)]
}

# common_dat = lapply(data_list, function(x) x[which(x$cpid %in% comcon_snps), ])
common_dat = map(data_list, ~ .x[which(.x$cpid %in% common_snps), ])

# Calculate some stats and plot:
median_se = unlist(lapply(common_dat, function(x) 1/(median(x$SE))))
c_value = unlist(lapply(common_dat, function(x) median(1/(sqrt(2 * x$MAF * (1-x$MAF))))))
sqrt_n = unlist(lapply(common_dat, function(x) sqrt(median(x$N))))

out_name = paste(out_dir, "all_SEvsN.png", sep = '')

png(out_name, width = 1000, height = 1000)
plot(x = c_value * median_se, y = sqrt_n, main = "SE vs. N", xlab = "C/median(SE)", ylab = "Sqrt(Nmax)", ylim = c(0, max(sqrt_n) * 1.1), xlim = c(0, max(c_value * median_se) * 1.1), pch = 20)
abline(0, 1, col = "red")
dev.off()

