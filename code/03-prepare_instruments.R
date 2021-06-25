rm(list=ls())

# Load libraries ---------------------------------------------------------------

library(magrittr)

# Specify paths ----------------------------------------------------------------

source("code/specify_paths.R")

# Load gwas information --------------------------------------------------------

risk_factors <- data.table::fread("data/risk_factors.csv", data.table = FALSE)

# Extract genome wide significant hits for each GWAS ---------------------------

df <- TwoSampleMR::extract_instruments(outcomes = unique(risk_factors$id),
                                       p1 = 5e-08,
                                       clump = TRUE,
                                       p2 = 5e-08,
                                       r2 = 0.001,
                                       kb = 10000)

colnames(df) <- gsub(".exposure","",colnames(df))

df[,c("mr_keep","pval_origin","id","data_source")] <- NULL

df$exposure <- gsub(".*id:","",df$exposure)
df$ncase <- NA
df$ncontrol <- NA

# Load outcome instruments -----------------------------------------------------

tmp <- data.table::fread("data/instruments_outcomes.csv")
df <- rbind(df, tmp)

# Make all effects positive ---------------------------------------------------

df$beta.original <- df$beta
df$eaf.original <- df$eaf
df$effect_allele.original <- df$effect_allele
df$other_allele.original <- df$other_allele
df[,c("beta","eaf","effect_allele","other_allele")] <- NULL

df <- df %>%
  dplyr::mutate(beta = ifelse(sign(beta.original)==-1, -1*beta.original, beta.original)) %>%
  dplyr::mutate(effect_allele = ifelse(sign(beta.original)==-1, other_allele.original, effect_allele.original)) %>%
  dplyr::mutate(other_allele = ifelse(sign(beta.original)==-1, effect_allele.original, other_allele.original)) %>%
  dplyr::mutate(eaf = ifelse(beta.original <= 0, 1-eaf.original, eaf.original))

df[,c("beta.original","eaf.original","effect_allele.original","other_allele.original")] <- NULL

# Save instruments -------------------------------------------------------------

data.table::fwrite(df,"data/instruments.csv")