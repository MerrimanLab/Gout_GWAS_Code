# Loading libraries.
library(fs)
library(vroom)
library(tidyverse)

# set colour palette
#devtools::install_github("G-Thomson/Manu")
library(Manu)

gwas_palette <- c("#7ACCD7", "#115896", "#7C6C65", "#4C4C53", "#BA2F00", "#B865A1", "#6AA086")
gwas_palette <- sort(gwas_palette)

# Setting path variable for working directory.
path1 <- path("/Volumes/archive/merrimanlab/copied/Nick/Major_GWAS/YY_Loci")

# Setting working directory.
setwd(path1)

# Loading urate GWAS summary stats column names.
urate_header <- read_delim(path(path1, "urate_header.txt"))

# Loading gout GWAS summary stats column names.
gout_header <- read_delim(path(path1, "gout_header.txt"))

# Loading in YY loci list urate effects.
yy_loci_urate <- read_delim(path(path1, "yy_snps_urate.txt"), col_names = F)
colnames(yy_loci_urate) <- colnames(urate_header)

# Loading in YY loci list gout effects.
yy_loci_gout <- read_delim(path(path1, "yy_snps_gout.txt"), col_names = F)
colnames(yy_loci_gout) <- colnames(gout_header)

# Loading in H4 loci list urate effects.
h4_loci_urate <- read_delim(path(path1, "h4_snps_urate.txt"), col_names = F)
colnames(h4_loci_urate) <- colnames(urate_header)

# Loading in H4 loci list gout effects.
h4_loci_gout <- read_delim(path(path1, "h4_snps_gout.txt"), col_names = F)
colnames(h4_loci_gout) <- colnames(gout_header)

# Extracting H4 urate vs gout effects for regression.
h4_loci_urate2 <- h4_loci_urate %>%
  select(CHR:POS, SNP, allele1:freq1, effect:P, N) %>%
  rename(EA_urate = allele1,
         AA_urate = allele2,
         EAF_urate = freq1,
         Beta_urate = effect,
         SE_urate = SE,
         P_urate = P,
         N_urate = N,
         BP = POS)

h4_loci_gout2 <- h4_loci_gout %>%
  select(CHR:SNP, minor:MAF, effect:P, N) %>%
  rename(EA_gout = minor,
         AA_gout = major,
         EAF_gout = MAF,
         Beta_gout = effect,
         SE_gout = SE,
         P_gout = P,
         N_gout = N)

h4_loci <- left_join(h4_loci_urate2, h4_loci_gout2)

# Ensuring that the variants are the same between GWAS using alleles and MAF.
h4_loci2 <- h4_loci %>%
  filter(EA_urate == EA_gout & AA_urate == AA_gout) %>%
  rename(EA = EA_urate,
         AA = AA_urate) %>%
  select(-EA_gout, -AA_gout)

tmp <- h4_loci2 %>%
  mutate(diff = EAF_gout - EAF_urate) %>%
  arrange(desc(abs(diff)))

# Flipping alleles to always be urate risk alleles.
h4_loci3 <- h4_loci2 %>%
  mutate(Beta_gout = case_when(Beta_urate < 0 ~ -1 * Beta_gout,
                               TRUE ~ Beta_gout),
         EAF_gout = case_when(Beta_urate < 0 ~ 1 - EAF_gout,
                              TRUE ~ EAF_gout),
         EAF_urate = case_when(Beta_urate < 0 ~ 1 - EAF_urate,
                               TRUE ~ EAF_urate),
         EA_tmp = EA,
         EA = case_when(Beta_urate < 0 ~ AA,
                        TRUE ~ EA),
         AA = case_when(Beta_urate < 0 ~ EA_tmp,
                        TRUE ~ AA),
         Beta_urate = case_when(Beta_urate < 0 ~ -1 * Beta_urate,
                                TRUE ~ Beta_urate)) %>%
  select(-EA_tmp)

# Modeling urate vs gout effects.
mod <- lm(Beta_gout ~ Beta_urate, data = h4_loci3)

# Extracting YY effects.
yy_loci_urate2 <- yy_loci_urate %>%
  select(CHR:POS, SNP, allele1:freq1, effect:P, N) %>%
  rename(EA_urate = allele1,
         AA_urate = allele2,
         EAF_urate = freq1,
         Beta_urate = effect,
         SE_urate = SE,
         P_urate = P,
         N_urate = N,
         BP = POS)

yy_loci_gout2 <- yy_loci_gout %>%
  select(CHR:SNP, minor:MAF, effect:P, N) %>%
  rename(EA_gout = minor,
         AA_gout = major,
         EAF_gout = MAF,
         Beta_gout = effect,
         SE_gout = SE,
         P_gout = P,
         N_gout = N)

yy_loci_effects <- left_join(yy_loci_urate2, yy_loci_gout2)

# Ensuring they are the same variants using alleles and EAF (removes 1 locus).
yy_loci_effects2 <- yy_loci_effects %>%
  filter(EA_urate == EA_gout & AA_urate == AA_gout) %>%
  rename(EA = EA_urate,
         AA = AA_urate) %>%
  select(-EA_gout, -AA_gout)

yy_loci_effects3 <- yy_loci_effects2 %>%
  mutate(diff = EAF_gout - EAF_urate) %>%
  filter(diff < 0.1)

# Flipping alleles to always be urate risk alleles.
yy_loci_effects4 <- yy_loci_effects3 %>%
  mutate(Beta_gout = case_when(Beta_urate < 0 ~ -1 * Beta_gout,
                               TRUE ~ Beta_gout),
         EAF_gout = case_when(Beta_urate < 0 ~ 1 - EAF_gout,
                              TRUE ~ EAF_gout),
         EAF_urate = case_when(Beta_urate < 0 ~ 1 - EAF_urate,
                               TRUE ~ EAF_urate),
         EA_tmp = EA,
         EA = case_when(Beta_urate < 0 ~ AA,
                        TRUE ~ EA),
         AA = case_when(Beta_urate < 0 ~ EA_tmp,
                        TRUE ~ AA),
         Beta_urate = case_when(Beta_urate < 0 ~ -1 * Beta_urate,
                                TRUE ~ Beta_urate)) %>%
  select(-EA_tmp)

# Removing low MAF variants (< 1%) - removes 3 variants.
yy_loci_effects5 <- yy_loci_effects4 %>%
  filter(EAF_urate > 0.01)

# Excluding variants with N < 1,000,000 for gout GWAS (removes 10 variants).
yy_loci_effects6 <- yy_loci_effects5 |>
  filter(N_gout > 1000000)

# Calculating expected gout effect for these 38 filtered YY loci.
yy_loci_effects7 <- yy_loci_effects6 %>%
  mutate(Beta_gout_exp = predict(mod, newdata = yy_loci_effects6))

# Making function to calculate power.
power_calc <- function(n, phi, b, f, threshold) {
  power <- pchisq(qchisq(threshold, df = 1, lower = F), df = 1, ncp = 2*f*(1-f)*n*phi*(1-phi)*b^2, lower = F)

  return(power)
}

# Calculating case:control ratio based on the European GWAS N cases vs N controls.
phi <- 100661/2106003

# Calculating total N.
n <- 100661 + 2106003

# Flipping effects to represent minor allele only.
yy_minor <- yy_loci_effects7 |>
  select(SNP, EAF_gout, Beta_gout_exp) |>
  mutate(MAF = case_when(EAF_gout > 0.5 ~ 1 - EAF_gout,
                         TRUE ~ EAF_gout))

# Extracting expected gout effects.
b <- yy_minor %>%
  pull(Beta_gout_exp)

# Extracting minor allele frequencies.
f <- yy_minor %>%
  pull(MAF)

# Calculating power for each variant at alpha of 0.01.
tmp <- c()
for(i in 1:length(b)) {
  tmp[i] <- power_calc(n, phi, b[i], f[i], 0.01)
}

# Extracting minimum beta where power > 99% for all MAFs at alpha = 0.01.
f_all <- seq(0.01, 0.49, by = 0.001)
b_all <- seq(0.001, 0.299, by = 0.001)
p_all <- c(5e-8, 1e-6, 1e-4, 0.01, 0.05)

min_b_all <- c()
for(k in 1:length(p_all)) {
  min_b <- c()
  for(j in 1:length(f_all)) {
    pwr_all <- c()
    for(i in 1:length(b_all)) {
      pwr_all[i] <- power_calc(n, phi, b_all[i], f_all[j], p_all[k])
    }

    if(max(pwr_all) < 0.99) {
      min_b[j] <- NA
    } else {
      a <- 1
      x <- 0
      while(x < 0.99) {
        a <- a + 1
        x <- pwr_all[a]
      }

      min_b[j] <- b_all[a]
    }
  }

  min_b_all <- c(min_b_all, min_b)
}

dat2 <- tibble("Effect size (logOR)" = min_b_all,
               "MAF" = rep(f_all, 5),
               "P" = rep(p_all, each = length(f_all)))
