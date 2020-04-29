rm(list=setdiff(ls(), keep))
graphics.off()

# Source functions -------------------------------------------------------------
 
source("code/fn-create_analysis_data.R")
source("code/fn-extract_feature.R")

# Load features data -----------------------------------------------------------

features <- readxl::read_xlsx("raw/feature_sources_ieugwas.xlsx",sheet = "All")
features <- features[!is.na(features$ieugwas),]

# Load type 2 diabetes GWAS ----------------------------------------------------

t2d <- data.table::fread("data/t2d_rsid.txt",
                         stringsAsFactors = FALSE,
                         data.table = FALSE)

t2d$samplesize <- NULL

# Load type 2 diabetes instrument ----------------------------------------------

t2d_ins <- data.table::fread("data/instrument-t2d.txt",
                             stringsAsFactors = FALSE,
                             data.table = FALSE)

# Perform analysis for each outcome --------------------------------------------

out = "cad"

# Load outcome GWAS ------------------------------------------------------------

outcome <- data.table::fread(paste0("data/",out,".txt"),
                             stringsAsFactors = FALSE,
                             data.table = FALSE)

# Create results data frame ----------------------------------------------------

# results <- data.frame(analysis = character(),
#                       direct = logical(),
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

# Perform analysis for each feature --------------------------------------------

for (i in 38:nrow(features)) {
  
  print(paste0("Row: ",i,"; Trait: ",features$trait[i]))
  
  ## Extract feature SNPs --------------------------------------------------------
  
  analysis_feature <- extract_feature(features$trait[i],features$ieugwas[i])
  
  ## Extract t2d SNPs ------------------------------------------------------------
  
  analysis_t2d <- t2d[t2d$SNP %in% analysis_feature$SNP,]
  
  ## Extract outcome SNPs --------------------------------------------------------
  
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
    results[nrow(results)+1,] <- c("MVMR",TRUE,features$trait[i],"type 2 diabetes",out,mvmr_out$coef[1,c(1,2,4)],length(unique(mvmr_in$SNP)))
    results[nrow(results)+1,] <- c("MVMR",TRUE,features$trait[i],features$trait[i],out,mvmr_out$coef[2,c(1,2,4)],length(unique(mvmr_in$SNP)))
    
    ## Perform univariate MR for feature -------------------------------------------
    
    uvmr_in <- MendelianRandomization::mr_input(bx = df[df$SNP %in% analysis_feature[analysis_feature$instrument==analysis_feature$trait,]$SNP,]$beta.feature,
                                                bxse = df[df$SNP %in% analysis_feature[analysis_feature$instrument==analysis_feature$trait,]$SNP,]$se.feature,
                                                by = df[df$SNP %in% analysis_feature[analysis_feature$instrument==analysis_feature$trait,]$SNP,]$beta.outcome,
                                                byse = df[df$SNP %in% analysis_feature[analysis_feature$instrument==analysis_feature$trait,]$SNP,]$se.outcome)
    uvmr_out <- MendelianRandomization::mr_ivw(uvmr_in)
    results[nrow(results)+1,] <- c("UVMR",NA,features$trait[i],features$trait[i],out,uvmr_out@Estimate,uvmr_out@StdError,uvmr_out@Pvalue,length(unique(analysis_feature[analysis_feature$instrument==analysis_feature$trait,]$SNP)))
    
    ## Perform univariate MR for type 2 diabetes -----------------------------------
    
    uvmr_in <- MendelianRandomization::mr_input(bx = df[df$SNP %in% t2d_ins$SNP,]$beta.t2d,
                                                bxse = df[df$SNP %in% t2d_ins$SNP,]$se.t2d,
                                                by = df[df$SNP %in% t2d_ins$SNP,]$beta.outcome,
                                                byse = df[df$SNP %in% t2d_ins$SNP,]$se.outcome)
    uvmr_out <- MendelianRandomization::mr_ivw(uvmr_in)
    results[nrow(results)+1,] <- c("UVMR",NA,features$trait[i],"type 2 diabetes",out,uvmr_out@Estimate,uvmr_out@StdError,uvmr_out@Pvalue,length(unique(t2d_ins$SNP)))
    
    ## Calculate indirect effects --------------------------------------------------
    
    ### Type 2 diabetes ------------------------------------------------------------
    
    direct_beta.t2d <- results[results$analysis=="MVMR" & results$feature==features$trait[i] & results$exposure=="type 2 diabetes",]$beta
    univar_beta.t2d <- results[results$analysis=="UVMR" & results$feature==features$trait[i] & results$exposure=="type 2 diabetes",]$beta
    
    direct_se.t2d <- results[results$analysis=="MVMR" & results$feature==features$trait[i] & results$exposure=="type 2 diabetes",]$se
    univar_se.t2d <- results[results$analysis=="UVMR" & results$feature==features$trait[i] & results$exposure=="type 2 diabetes",]$se
    
    indirect_beta.t2d <- as.numeric(univar_beta.t2d) - as.numeric(direct_beta.t2d)
    indirect_se.t2d <- sqrt((as.numeric(univar_se.t2d))^2 + as.numeric((direct_se.t2d))^2)
    results[nrow(results)+1,] <- c("MVMR",FALSE,features$trait[i],"type 2 diabetes",out,indirect_beta.t2d,indirect_se.t2d,NA,NA)
    
    ### Feature --------------------------------------------------------------------
    
    direct_beta.feature <- results[results$analysis=="MVMR" & results$feature==features$trait[i] & results$exposure==features$trait[i],]$beta
    univar_beta.feature <- results[results$analysis=="UVMR" & results$feature==features$trait[i] & results$exposure==features$trait[i],]$beta
    
    direct_se.feature <- results[results$analysis=="MVMR" & results$feature==features$trait[i] & results$exposure==features$trait[i],]$se
    univar_se.feature <- results[results$analysis=="UVMR" & results$feature==features$trait[i] & results$exposure==features$trait[i],]$se
    
    indirect_beta.feature  <- as.numeric(univar_beta.feature) - as.numeric(direct_beta.feature)
    indirect_se.feature <- sqrt((as.numeric(univar_se.feature))^2 + as.numeric((direct_se.feature))^2)
    results[nrow(results)+1,] <- c("MVMR",FALSE,features$trait[i],features$trait[i],out,indirect_beta.feature,indirect_se.feature,NA,NA)
    
  } else {
    
    results[nrow(results)+1,] <- c("MVMR",TRUE,features$trait[i],"type 2 diabetes",out,rep(NA,4))
    results[nrow(results)+1,] <- c("MVMR",TRUE,features$trait[i],features$trait[i],out,rep(NA,4))
    results[nrow(results)+1,] <- c("UVMR",NA,features$trait[i],"type 2 diabetes",out,rep(NA,3),length(unique(t2d_ins$SNP)))
    results[nrow(results)+1,] <- c("UVMR",NA,features$trait[i],features$trait[i],out,rep(NA,3),length(unique(analysis_feature[analysis_feature$instrument==analysis_feature$trait,]$SNP)))
    results[nrow(results)+1,] <- c("MVMR",FALSE,features$trait[i],"type 2 diabetes",out,rep(NA,4))
    results[nrow(results)+1,] <- c("MVMR",FALSE,features$trait[i],features$trait[i],out,rep(NA,4))
    
  }
  
  data.table::fwrite(results,"output/mvmr_results.csv",row.names = FALSE)
  
}
