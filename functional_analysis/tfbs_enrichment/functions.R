# Functions for doing TFBS analysis

# Function to align the meQTL effect size to the gout risk allele
align_effect <- function (data) {
  ind = which(data$effect.gwas < 0)
  data$effect.gwas[ind] = -data$effect.gwas[ind]
  data$MAF.gwas[ind] = 1 - data$MAF.gwas[ind]
  data$effect.mqtl[ind] = -data$effect.mqtl[ind]
  data$MAF.mqtl[ind] = 1 - data$MAF.mqtl[ind]
  tmp = data$minor[ind]
  data$minor[ind] = data$major[ind]
  data$major[ind] = tmp
  return(data)
}

# Function to rbind all of the common TFBS data from a list
simplify_tf_dat <- function (data_list, index) {
  tf_list = unique(index$TF)
  res = map(tf_list, ~ map_dfr(.x, ~ data_list[which(index$TF == .x)] %>% map(., ~ .x[,1:3])) %>% distinct)
  names(res) = tf_list
  return(res)
}

# Function to load TFBS data of a specific cell type
pull_tfbs <- function (chip_info, group, ignore.group = F) {
  if (ignore.group) {
    cell_tfbs = chip_info
  } else {
    cell_tfbs = chip_info %>% filter(Group == group)
  }
  cell_dat = cell_tfbs %>% pull(label) %>% paste('data/RELI_data/ChIP-seq/', ., sep = '') %>% map(., ~ read.table(.x, header = F, sep = '\t', stringsAsFactors = F)) %>% map(., ~ .x %>% mutate(V1 = gsub('chr', '', V1, ignore.case = T)))
  cell_dat = simplify_tf_dat(cell_dat, cell_tfbs)
  return(cell_dat)
}

# Function to return a vector of TRUE/FALSE of whether a site overlapped with
# a given TFBS data
check_overlap <- function (query, tf_dat, offset = 50) {
  res = c()
  for (i in 1:nrow(query)) {
    chr = tf_dat$V1
    start = tf_dat$V2 - offset
    end = tf_dat$V3 + offset
    tmp = any(chr == query$CHR[i] & start <= query$BP[i] & end >= query$BP[i])
    res = c(res, tmp)
  }
  return(res)
}

# Function to fill the TRUE value with effect size, while FALSE value with 0
fill_effect <- function (effect_dat, tf_dat) {
  for (i in 1:ncol(tf_dat)) {
    ind1 = which(tf_dat[, i])
    ind2 = which(!tf_dat[, i])
    tf_dat[ind1, i] = effect_dat$effect.mqtl[ind1]
    tf_dat[ind2, i] = 0
  }
  return(tf_dat)
}

# Function to generate heatmap for a specific cell-type TFBS
run_heatmap <- function (data, filename = 'tmp.png') {
  png(filename, height = 1000, width = 1000)
  heatmap.2(data, col = 'bluered', trace = 'none', scale = 'none')
  dev.off()
}

# Function to generate a heatmap-ready data, given a list of TF meta data.
# It will output a list where each list element is a TF by CpG matrix
compile_data <- function (chip_info, cpg_dat, offset = 50, ignore.group = F) {
  cpg_chrpos = cpg_dat %>% mutate(CHR = CHR.mqtl, BP = BP.mqtl)
  snp_cpg_names = paste(cpg_dat$SNP, cpg_dat$cpg, sep = '_')
  if (ignore.group) {
    cell_types = 'ALL'
    tf_list = pull_tfbs(chip_info, NULL, ignore.group = T)
    res = map_dfc(tf_list, ~ check_overlap(cpg_chrpos, .x, offset)) %>% as.data.frame
    res = fill_effect(cpg_chrpos, res) %>% t() %>% as.matrix
    colnames(res) = snp_cpg_names
  } else {
    cell_types = unique(chip_info$Group)
    tf_list = map(cell_types, ~ pull_tfbs(chip_info, .x))
    res = map(tf_list, ~ map_dfc(.x, ~ check_overlap(cpg_chrpos, .x, offset)) %>% as.data.frame)
    res = map(res, ~ fill_effect(cpg_chrpos, .x))
    res = map(res, ~ as.matrix(t(.x)))
    res = map(res, .f = function(x) {colnames(x) = snp_cpg_names; return(x)})
  }
  names(res) = cell_types
  return(res)
}

# Function to generate contingency tables for chi-square test
#
# | A   | B   || X
# | --- | --- || --
# | C   | D   || -
# | === | === || ==
# | Y   | -   || E
#
# count1 = n(TF binding in set1) (A)
# count2 = n(TF binding in set2) (X)
# set_n = total number of values in set1 (Y)
# total = total number of values used in the table (E)
make_table <- function (count1, count2, set_n = 402, total = 232477) {
  val1 = count1
  val2 = count2 - val1
  val3 = set_n - count1
  val4 = total - sum(c(val1, val2, val3))
  tab = matrix(c(val1, val2, val3, val4), c(2, 2), byrow = T)
  return(tab)
}

# Function to run Fisher's exact or chi-square test, given a 2 by 2 contingency
# table
# return_p = only return p-value
test_table <- function (data, test = 'fisher', return_p = T, return_sig = F) {
  if (test == 'chi') {
    if (any(data == 0)) {
      stop('There is a 0 in one of the cells in the contingency table.')
    }
    # Simulate the p-value if there is a cell with < 20 counts. Otherwise run
    # default chi-square test
    test = ifelse(any(data < 20), chisq.test(data, simulate.p.value = T, B = 1e6), chisq.test(data))
    p = ifelse(any(data < 20), chisq.test(data, simulate.p.value = T, B = 1e6)$p.value, chisq.test(data)$p.value)
  } else if (test == 'fisher') {
    test = fisher.test(data)
    p = test$p.value
  }
  if (return_p) {
    return(p)
  }
  return(test)
}

