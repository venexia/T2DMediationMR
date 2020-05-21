rm(list=ls())
graphics.off()

# Load packages ----------------------------------------------------------------

library(tidyverse)

# Specify paths ----------------------------------------------------------------

source("code/specify_paths.R")

# Source functions -------------------------------------------------------------

source("code/fn-rsid.R")
source("code/fn-extract_feature_ieugwas.R")
source("code/fn-extract_feature_manual.R")
source("code/fn-create_analysis_data.R")
source("code/fn-perform_uvmr.R")

# Load features data -----------------------------------------------------------

features <- data.table::fread("data/features_sources.csv",
                              stringsAsFactors = FALSE,
                              data.table = FALSE)

features <- features[!(features$trait %in% c("testosterone_f","testosterone_m")),]

# Load type 2 diabetes instrument ----------------------------------------------

t2d_ins <- data.table::fread("data/instrument-t2d.txt",
                             stringsAsFactors = FALSE,
                             data.table = FALSE)

# Load outcome GWAS ------------------------------------------------------------

pad <- data.table::fread("data/pad.txt",
                         stringsAsFactors = FALSE,
                         data.table = FALSE)

cad <- data.table::fread("data/cad.txt",
                         stringsAsFactors = FALSE,
                         data.table = FALSE)

# Create results data frame ----------------------------------------------------

results <- data.frame(exposure = character(),
                      outcome = character(),
                      estimate = numeric(),
                      se = numeric(),
                      pval = numeric(),
                      snps = character(),
                      stringsAsFactors = FALSE)

# Perform analysis for each feature-outcome combination ------------------------

for (i in c("cad","pad")) {
  
  out <- get(i)
  
  ## Create analysis dataset --------------------------------------------------------
  
  df <- t2d_ins[,c("SNP","effect_allele","other_allele","beta","se")]
  colnames(df) <- c("SNP","effect_allele.exposure","other_allele.exposure","beta.exposure","se.exposure")
  
  df <- merge(df,out[,c("SNP","effect_allele","other_allele","beta","se")],all.x = TRUE,by = "SNP")
  colnames(df) <- c("SNP","effect_allele","other_allele","beta.exposure","se.exposure",
                    "effect_allele.outcome","other_allele.outcome","beta.outcome.orig","se.outcome")
  
  ## Harmonise direction of effect -----------------------------------------------
  
  df$beta.outcome <- NA
  df$beta.outcome <- ifelse(df$effect_allele==df$effect_allele.outcome & df$other_allele==df$other_allele.outcome,df$beta.outcome.orig,df$beta.outcome)
  df$beta.outcome <- ifelse(df$effect_allele==df$other_allele.outcome & df$other_allele==df$effect_allele.outcome,-1*df$beta.outcome.orig,df$beta.outcome)
  
  df <- na.omit(unique(df))
  
  df <- df[,c("SNP","effect_allele","other_allele","beta.exposure","se.exposure","beta.outcome","se.outcome")]
  
  # Prepare MR dataset
  
  uvmr_in <- MendelianRandomization::mr_input(bx = df$beta.exposure,
                                              bxse = df$se.exposure,
                                              by = df$beta.outcome,
                                              byse = df$se.outcome)
  
  # Perform MR
  
  uvmr_out <- MendelianRandomization::mr_ivw(uvmr_in)
  
  # Save results
  
  results[nrow(results)+1,] <- c("type 2 diabetes",i,uvmr_out@Estimate,uvmr_out@StdError,uvmr_out@Pvalue,paste0(df$SNP,collapse = ";"))
  
}

results$lci <- as.numeric(results$estimate) - qnorm(0.975)*as.numeric(results$se)
results$uci <- as.numeric(results$estimate) + qnorm(0.975)*as.numeric(results$se)
results$or <- exp(as.numeric(results$estimate))
results$lci_or <- exp(as.numeric(results$lci))
results$uci_or <- exp(as.numeric(results$uci))

data.table::fwrite(results,"output/uvmr_t2d_outcomes.csv",row.names = FALSE)