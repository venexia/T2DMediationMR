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

# Create results data frame ----------------------------------------------------

# results <- data.frame(effect = character(),
#                       feature = character(),
#                       exposure = character(),
#                       outcome = character(),
#                       beta = numeric(),
#                       se = numeric(),
#                       pval = numeric(),
#                       nsnp = numeric(),
#                       stringsAsFactors = FALSE)

results <- data.table::fread("output/mvmr_results.csv",
                             data.table = FALSE,
                             stringsAsFactors = FALSE)

# Perform analysis for each outcome --------------------------------------------

for (out in c("cad")) { # "pad"
  
  # Load outcome GWAS ------------------------------------------------------------
  
  outcome <- data.table::fread(paste0("data/",out,".txt"),
                               stringsAsFactors = FALSE,
                               data.table = FALSE)

  # Perform analysis for each feature --------------------------------------------
  
  for (i in 1:nrow(features)) {
    
    print(paste0("Outcome: ",out,"; Trait: ",features$trait[i]))
    
    ## Extract feature SNPs --------------------------------------------------------
    
    if (features$ieugwas[i]!="") {
      analysis_feature <- extract_feature_ieugwas(features$trait[i],features$ieugwas[i])
    } else {
      analysis_feature <- extract_feature_manual(features$trait[i])
    }
    
    ## Extract snplist from T2D GWAS ---------------------------------------------
    
    analysis_t2d <- t2d[t2d$SNP %in% analysis_feature$SNP,]
    
    ## Extract snplist from outcome GWAS -----------------------------------------
    
    analysis_outcome <- outcome[outcome$SNP %in% analysis_feature$SNP,]
    
    ## Create clumped analysis data ------------------------------------------------
    
    df <- create_analysis_data(data_t2d = analysis_t2d,
                               data_feature = analysis_feature,
                               data_outcome = analysis_outcome)
    
    if (nrow(df[df$SNP %in% analysis_feature[analysis_feature$instrument==analysis_feature$trait,]$SNP,])>0) {
      
      ## Perform MVMR ----------------------------------------------------------------
      
      mvmr_in <- MVMR::format_mvmr(BXGs = df[,c("beta.t2d","beta.feature")],
                                   BYG = df[,c("beta.outcome")],
                                   seBXGs = df[,c("se.t2d","se.feature")],
                                   seBYG = df[,c("se.outcome")],
                                   RSID = df[,c("SNP")])
      mvmr_out <- MVMR::mvmr(r_input = mvmr_in, gencov = 0, weights = 1)
      results[nrow(results)+1,] <- c("direct",features$trait[i],"type 2 diabetes",out,mvmr_out$coef[1,c(1,2,4)],length(unique(mvmr_in$SNP)))
      results[nrow(results)+1,] <- c("direct",features$trait[i],features$trait[i],out,mvmr_out$coef[2,c(1,2,4)],length(unique(mvmr_in$SNP)))
      
      ## Perform univariate MR for feature -------------------------------------------
      
      uvmr_in <- MendelianRandomization::mr_input(bx = df[df$SNP %in% analysis_feature[analysis_feature$instrument==analysis_feature$trait,]$SNP,]$beta.feature,
                                                  bxse = df[df$SNP %in% analysis_feature[analysis_feature$instrument==analysis_feature$trait,]$SNP,]$se.feature,
                                                  by = df[df$SNP %in% analysis_feature[analysis_feature$instrument==analysis_feature$trait,]$SNP,]$beta.outcome,
                                                  byse = df[df$SNP %in% analysis_feature[analysis_feature$instrument==analysis_feature$trait,]$SNP,]$se.outcome)
      uvmr_out <- MendelianRandomization::mr_ivw(uvmr_in)
      results[nrow(results)+1,] <- c("total",features$trait[i],features$trait[i],out,uvmr_out@Estimate,uvmr_out@StdError,uvmr_out@Pvalue,length(unique(analysis_feature[analysis_feature$instrument==analysis_feature$trait,]$SNP)))
      
      ## Perform univariate MR for type 2 diabetes -----------------------------------
      
      uvmr_in <- MendelianRandomization::mr_input(bx = df[df$SNP %in% t2d_ins$SNP,]$beta.t2d,
                                                  bxse = df[df$SNP %in% t2d_ins$SNP,]$se.t2d,
                                                  by = df[df$SNP %in% t2d_ins$SNP,]$beta.outcome,
                                                  byse = df[df$SNP %in% t2d_ins$SNP,]$se.outcome)
      uvmr_out <- MendelianRandomization::mr_ivw(uvmr_in)
      results[nrow(results)+1,] <- c("total",features$trait[i],"type 2 diabetes",out,uvmr_out@Estimate,uvmr_out@StdError,uvmr_out@Pvalue,length(unique(t2d_ins$SNP)))
      
      ## Calculate indirect effects --------------------------------------------------
      
      ### Type 2 diabetes ------------------------------------------------------------
      
      direct_beta.t2d <- results[results$effect=="direct" & results$feature==features$trait[i] & results$exposure=="type 2 diabetes",]$beta
      total_beta.t2d <- results[results$effect=="total" & results$feature==features$trait[i] & results$exposure=="type 2 diabetes",]$beta
      
      direct_se.t2d <- results[results$effect=="direct" & results$feature==features$trait[i] & results$exposure=="type 2 diabetes",]$se
      total_se.t2d <- results[results$effect=="total" & results$feature==features$trait[i] & results$exposure=="type 2 diabetes",]$se
      
      indirect_beta.t2d <- as.numeric(total_beta.t2d) - as.numeric(direct_beta.t2d)
      indirect_se.t2d <- sqrt((as.numeric(total_se.t2d))^2 + as.numeric((direct_se.t2d))^2)
      results[nrow(results)+1,] <- c("indirect",features$trait[i],"type 2 diabetes",out,indirect_beta.t2d,indirect_se.t2d,NA,NA)
      
      ### Feature --------------------------------------------------------------------
      
      direct_beta.feature <- results[results$effect=="direct" & results$feature==features$trait[i] & results$exposure==features$trait[i],]$beta
      total_beta.feature <- results[results$effect=="total" & results$feature==features$trait[i] & results$exposure==features$trait[i],]$beta
      
      direct_se.feature <- results[results$effect=="direct" & results$feature==features$trait[i] & results$exposure==features$trait[i],]$se
      total_se.feature <- results[results$effect=="total" & results$feature==features$trait[i] & results$exposure==features$trait[i],]$se
      
      indirect_beta.feature  <- as.numeric(total_beta.feature) - as.numeric(direct_beta.feature)
      indirect_se.feature <- sqrt((as.numeric(total_se.feature))^2 + as.numeric((direct_se.feature))^2)
      results[nrow(results)+1,] <- c("indirect",features$trait[i],features$trait[i],out,indirect_beta.feature,indirect_se.feature,NA,NA)
      
    } else {
      
      results[nrow(results)+1,] <- c("direct",features$trait[i],"type 2 diabetes",out,rep(NA,4))
      results[nrow(results)+1,] <- c("direct",features$trait[i],features$trait[i],out,rep(NA,4))
      results[nrow(results)+1,] <- c("total",features$trait[i],"type 2 diabetes",out,rep(NA,3),length(unique(t2d_ins$SNP)))
      results[nrow(results)+1,] <- c("total",features$trait[i],features$trait[i],out,rep(NA,3),length(unique(analysis_feature[analysis_feature$instrument==analysis_feature$trait,]$SNP)))
      results[nrow(results)+1,] <- c("indirect",features$trait[i],"type 2 diabetes",out,rep(NA,4))
      results[nrow(results)+1,] <- c("indirect",features$trait[i],features$trait[i],out,rep(NA,4))
      
    }
    
    data.table::fwrite(results,"output/mvmr_results.csv",row.names = FALSE)
    
  }
}