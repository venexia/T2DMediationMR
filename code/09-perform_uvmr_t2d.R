rm(list=ls())
graphics.off()

# Load packages ----------------------------------------------------------------

library(tidyverse)

# Specify paths ----------------------------------------------------------------

source("code/specify_paths.R")

# Source functions -------------------------------------------------------------

source("code/fn-extract_feature_ieugwas.R")
source("code/fn-extract_feature_manual.R")

# Load features data -----------------------------------------------------------

features <- data.table::fread("data/features_sources.csv",
                              stringsAsFactors = FALSE,
                              data.table = FALSE)

features <- features[!(features$trait %in% c("testosterone_f","testosterone_m")),]

# Load type 2 diabetes GWAS ----------------------------------------------------

t2d <- data.table::fread("data/t2d_rsid.txt",
                         stringsAsFactors = FALSE,
                         data.table = FALSE)


# Load type 2 diabetes instrument ----------------------------------------------

t2d_ins <- data.table::fread("data/instrument-t2d.txt",
                             stringsAsFactors = FALSE,
                             data.table = FALSE)

# Create results data frame ----------------------------------------------------

# results <- data.frame(exposure = character(),
#                       outcome = character(),
#                       estimate = numeric(),
#                       se = numeric(),
#                       pval = numeric(),
#                       snps = character(),
#                       stringsAsFactors = FALSE)
# 
# data.table::fwrite(results,"output/uvmr_t2d.csv",row.names = FALSE)

results <- data.table::fread("output/uvmr_t2d.csv",
                             stringsAsFactors = FALSE,
                             data.table = FALSE)

# Perform analysis for each feature-outcome combination ------------------------

for (i in 67:nrow(features)) {
  
  print(paste0("i: ",i))
  
  ## Extract feature SNPs --------------------------------------------------------
  
  if (features[i,]$ieugwas!="") {
    analysis_feature <- extract_feature_ieugwas(trait = features[i,]$trait,
                                                id = features[i,]$ieugwas)
  } else {
    analysis_feature <- extract_feature_manual(features[i,]$trait)
  }

  analysis_feature <- analysis_feature[analysis_feature$instrument!="type 2 diabetes",]
  
  if (nrow(analysis_feature)>0) {
    
  ## Clump feature SNP --------------------------------------------------------
  
  clumped <- analysis_feature[,c("SNP","pval")]
  colnames(clumped) <- c("rsid","pval")
  clumped <- ieugwasr::ld_clump(dat = clumped, clump_kb = 10000, clump_r2 = 0.001)

  ## Create analysis dataset --------------------------------------------------------
  
  df <- analysis_feature[analysis_feature$SNP %in% clumped$rsid,c("SNP","effect_allele","other_allele","beta","se")]
  colnames(df) <- c("SNP","effect_allele.exposure","other_allele.exposure","beta.exposure","se.exposure")
  
  df <- merge(df,t2d[,c("SNP","effect_allele","other_allele","beta","se")],all.x = TRUE,by = "SNP")
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
  
  results[nrow(results)+1,] <- c(features[i,]$trait,"type 2 diabetes",uvmr_out@Estimate,uvmr_out@StdError,uvmr_out@Pvalue,paste0(df$SNP,collapse = ";"))
  
  data.table::fwrite(results,"output/uvmr_t2d.csv",row.names = FALSE)
  
  } else {
    
    # Save results
    
    results[nrow(results)+1,] <- c(features[i,]$trait,"type 2 diabetes",rep(NA,3),"")
    
    data.table::fwrite(results,"output/uvmr_t2d.csv",row.names = FALSE)
    
  }
  
}

data.table::fwrite(results,"output/uvmr_t2d.csv",row.names = FALSE)
