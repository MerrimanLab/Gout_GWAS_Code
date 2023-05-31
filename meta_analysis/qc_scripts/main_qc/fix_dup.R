################################################################################
# Rscript to fix any remaining duplicate variants, after the datasets have gone
# through per_study_qc.R
################################################################################
library(dplyr)
library(vroom)
library(furrr)

args = commandArgs(trailingOnly = T)

data_list = list()

dat_names = gsub('\\..*', '', args)
dat_names = gsub('.*/', '', dat_names)

for (i in 1:length(args)) {
	data_list[[i]] = vroom(args[i])
}

names(data_list) = dat_names

reorder = order(unlist(lapply(data_list, nrow)), decreasing = T)

data_list = data_list[reorder]
dat_names = dat_names[reorder]

# Set up multi-process env for furrr:
plan(strategy = 'multicore', workers = length(data_list))
options(future.globals.maxSize= 50 * 1024 ^ 3 )

# Remove any "easy" duplicate (exact same allelic value):
message('Removing easy duplicates...')
# data_list = lapply(data_list, function(x) x %>% distinct(cpid, allelic_value, .keep_all = T))
data_list = future_map(data_list, ~ .x %>% distinct(cpid, allelic_value, .keep_all = T))

# Get duplicate cpid from each data set:
message('Identifying duplicated cpid...')
# dup_var = lapply(data_list, function(x) x$cpid[which(duplicated(x$cpid))])
dup_var = future_map(data_list, ~ .x$cpid[which(duplicated(.x$cpid))])

# Add CHR/POS/allelic value id to all data sets:
data_list = future_map(data_list, ~ .x %>% mutate(cpavid = paste(cpid, allelic_value, sep = '_')))

# Run a loop to pull allelic values for the variant in all data set and find
# the most concordant variant (i.e. find the most frequent allelic value across
# the data sets for that variant)
#
# NOTE: if the allelic value counts are tied (i.e. no max), this indicates that
# there is a mix of alleles across the different data sets for that position.
# For example, dat1 = "A/G", dat2 = "A/G" and "C/G" (i.e. duplicated chr/pos
# with different allele in this data set), and dat3 = "C/G". Inconsistent
# alleles will be removed in the next QC step, since the alleles don't match
# across the data sets. Therefore, just take the first variant and proceed.
message('Checking allele consistency across datasets...')
for (i in 1:length(dup_var)) {
	message('\t', i, ' out of ', length(dup_var))
	tmp_dup = dup_var[[i]]
	if (length(tmp_dup) == 0) next
	tmp_data_list = lapply(data_list, function(x) x[which(x$cpid %in% tmp_dup),])
	rm_list = c()
	for (j in 1:length(tmp_dup)) {
		tmp_dat = tmp_data_list[[i]]
		tmp_av = lapply(tmp_data_list, function(x) x$allelic_value[x$cpid == tmp_dup[j]]) %>% unlist %>% table
		ind_av = which(tmp_av == max(tmp_av))
		av = names(tmp_av)[ind_av[1]]
		rm_list = c(rm_list, tmp_dat$cpavid[(tmp_dat$cpid == tmp_dup[j] & tmp_dat$allelic_value != av)])
	}
	tmp_dat = data_list[[i]]
	tmp_dat = tmp_dat[which(!(tmp_dat$cpavid %in% rm_list)),]
	data_list[[i]] = tmp_dat
}

dat_name = gsub('\\..*', '', args[reorder])
file_name = paste(dat_name, "nodup.tsv", sep = '_')

# Save the non-duplicated data:
for (i in 1:length(data_list)) {
	message('Writing ', names(data_list[i]), ' as ', file_name[i])
	vroom_write(data_list[[i]], file_name[i])
}
