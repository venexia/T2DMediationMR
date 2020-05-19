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

# Load type 2 diabetes GWAS ----------------------------------------------------

t2d <- data.table::fread("data/t2d_rsid.txt",
                         stringsAsFactors = FALSE,
                         data.table = FALSE)

t2d$samplesize <- NULL

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

feature_outcome <- data.frame(feature = rep(unique(features$trait), times = 2),
                              outcome = rep(c("cad","pad"), each = length(unique(features$trait))),
                              stringsAsFactors = FALSE)

# results <- data.frame(analysis = character(),
#                       exposure = character(),
#                       outcome = character(),
#                       effect = character(),
#                       estimate = numeric(),
#                       se = numeric(),
#                       pval = numeric(),
#                       snps = character(),
#                       Q_strength = character(),
#                       Q_valid = character(),
#                       stringsAsFactors = FALSE)
# 
# data.table::fwrite(results,"output/mvmr_results.csv",row.names = FALSE)

results <- data.table::fread("output/mvmr_results.csv",
                             stringsAsFactors = FALSE,
                             data.table = FALSE)

# Perform analysis for each feature-outcome combination ------------------------

for (i in 59:nrow(feature_outcome)) {
  
  print(paste0("i: ",i,"; Feature: ",feature_outcome$feature[i],"; Outcome: ",feature_outcome$outcome[i]))
  
  analysis <- paste0(feature_outcome$feature[i],"_",feature_outcome$outcome[i])
  
  ## Extract feature SNPs --------------------------------------------------------
  
  if (features[features$trait==feature_outcome$feature[i],]$ieugwas!="") {
    analysis_feature <- extract_feature_ieugwas(features[features$trait==feature_outcome$feature[i],]$trait,features[features$trait==feature_outcome$feature[i],]$ieugwas)
  } else {
    analysis_feature <- extract_feature_manual(features[features$trait==feature_outcome$feature[i],]$trait)
  }
  
  ## Extract snplist from T2D GWAS ---------------------------------------------
  
  analysis_t2d <- t2d[t2d$SNP %in% analysis_feature$SNP,]
  
  ## Extract snplist from outcome GWAS -----------------------------------------
  
  outcome <- get(feature_outcome$outcome[i])
  analysis_outcome <- outcome[outcome$SNP %in% analysis_feature$SNP,]
  
  ## Create clumped analysis data ------------------------------------------------
  
  df <- create_analysis_data(data_t2d = analysis_t2d,
                             data_feature = analysis_feature,
                             data_outcome = analysis_outcome)
  
  instrument_t2d <- intersect(t2d_ins$SNP,df$SNP)
  instrument_feature <- intersect(analysis_feature[analysis_feature$instrument==feature_outcome$feature[i],]$SNP,df$SNP)
  
  if (nrow(df[df$SNP %in% analysis_feature[analysis_feature$instrument==analysis_feature$trait,]$SNP,])>0) {
    
    ## Perform MVMR ----------------------------------------------------------------
    
    mvmr_in <- MVMR::format_mvmr(BXGs = df[,c("beta.t2d","beta.feature")],
                                 BYG = df[,c("beta.outcome")],
                                 seBXGs = df[,c("se.t2d","se.feature")],
                                 seBYG = df[,c("se.outcome")],
                                 RSID = df[,c("SNP")])
    mvmr_out <- MVMR::mvmr(r_input = mvmr_in, gencov = 0, weights = 1)
    
    results[nrow(results)+1,] <- c(analysis,"t2d",feature_outcome$outcome[i],"direct",mvmr_out$coef[1,c(1,2,4)],paste0(mvmr_in$SNP,collapse = ";"),mvmr_out$Q_strength[1],mvmr_out$Q_valid)
    results[nrow(results)+1,] <- c(analysis,feature_outcome$feature[i],feature_outcome$outcome[i],"direct",mvmr_out$coef[2,c(1,2,4)],paste0(mvmr_in$SNP,collapse = ";"),mvmr_out$Q_strength[1],mvmr_out$Q_valid)
    
    ## Perform univariate MR ---------------------------------------------------
    
    results[nrow(results)+1,] <- perform_uvmr(analysis,feature_outcome$feature[i],"outcome",feature_outcome$outcome[i])                 
    results[nrow(results)+1,] <- perform_uvmr(analysis,"t2d","outcome",feature_outcome$outcome[i])
    results[nrow(results)+1,] <- perform_uvmr(analysis,feature_outcome$feature[i],"t2d","t2d")
    results[nrow(results)+1,] <- perform_uvmr(analysis,"t2d",feature_outcome$feature[i],feature_outcome$feature[i])
    
  } else {
    
    results[nrow(results)+1,] <- c(analysis,"t2d",feature_outcome$outcome[i],"direct",rep(NA,3),"",rep(NA,2))
    results[nrow(results)+1,] <- c(analysis,feature_outcome$feature[i],feature_outcome$outcome[i],"direct",rep(NA,3),"",rep(NA,2))
    results[nrow(results)+1,] <- c(analysis,feature_outcome$feature[i],feature_outcome$outcome[i],"total",rep(NA,3),"",rep(NA,2))                
    results[nrow(results)+1,] <- c(analysis,"t2d",feature_outcome$outcome[i],"total",rep(NA,3),"",rep(NA,2))     
    results[nrow(results)+1,] <- c(analysis,feature_outcome$feature[i],"t2d","total",rep(NA,3),"",rep(NA,2))     
    results[nrow(results)+1,] <- c(analysis,"t2d",feature_outcome$feature[i],"total",rep(NA,3),"",rep(NA,2))     
    
  }
  
  data.table::fwrite(results,"output/mvmr_results.csv",row.names = FALSE)
  
}