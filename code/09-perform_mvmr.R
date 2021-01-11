rm(list=ls())
graphics.off()

# Specify paths ----------------------------------------------------------------

source("code/specify_paths.R", echo = TRUE)

# Source functions -------------------------------------------------------------

source("code/fn-extract_snps.R", echo = TRUE)
source("code/fn-harmonise_snps.R", echo = TRUE)

# Specify results data frame ---------------------------------------------------

# results <- data.frame(analysis = character(),
#                       exposure = character(),
#                       outcome = character(),
#                       effect = character(),
#                       estimate = numeric(),
#                       se = numeric(),
#                       pval = numeric(),
#                       snps = character(),
#                       Q_strength = numeric(),
#                       Q_valid = numeric(),
#                       condF = numeric(),
#                       stringsAsFactors = FALSE)

results <- data.table::fread("output/mvmr_results.csv", data.table = FALSE)

# Specify type 2 diabetes instrument to use ------------------------------------

t2d_instrument <- "t2d"

# Load source data info --------------------------------------------------------

gwas <- data.table::fread("raw/gwas.csv",
                          stringsAsFactors = FALSE,
                          data.table = FALSE)

features <- gwas[gwas$feature==TRUE,]$trait

for (j in c("pad")) {

  # Load instruments -------------------------------------------------------------
  
  instruments <- data.table::fread("data/instruments_all.txt",
                                   data.table = FALSE)
  
  # instruments <- data.table::fread(paste0("data/instruments_",j,".txt"),
  #                                  data.table = FALSE)
  
  instruments <- instruments[instruments$exposure %in% c(features,t2d_instrument),]
  
  for (i in unique(instruments[instruments$exposure!=t2d_instrument,]$exposure)[1:45]) {
    
    analysis <- paste(i, j, sep = "_")
    print(paste0("Analysis: ",analysis))
  
    # List MVMR instrument SNPs ----------------------------------------------
    
    snplist <- unique(instruments[instruments$exposure %in% c(i,t2d_instrument),
                           c("SNP","chrpos")])
    
    snplist$instrument <- ""
    snplist$instrument <- ifelse(snplist$SNP %in% instruments[instruments$exposure==i,]$SNP & snplist$SNP %in% instruments[instruments$exposure==t2d_instrument,]$SNP,"both",snplist$instrument)
    snplist$instrument <- ifelse(snplist$SNP %in% instruments[instruments$exposure==i,]$SNP & !(snplist$SNP %in% instruments[instruments$exposure==t2d_instrument,]$SNP),"feature",snplist$instrument)
    snplist$instrument <- ifelse(!(snplist$SNP %in% instruments[instruments$exposure==i,]$SNP) & snplist$SNP %in% instruments[instruments$exposure==t2d_instrument,]$SNP,"t2d",snplist$instrument)

    # Extract MVMR instrument SNPs from GWAS ---------------------------------
    
    t2d <- extract_snps(gwas = "t2d",
                        ieugwas = "",
                        snplist = snplist,
                        label = "t2d")
    
    feature <- extract_snps(gwas = i,
                            ieugwas = gwas[gwas$trait==i,]$ieugwas,
                            snplist = snplist,
                            label = "feature")
    
    outcome <- extract_snps(gwas = j,
                            ieugwas = gwas[gwas$trait==j,]$ieugwas,
                            snplist = snplist,
                            label = "outcome")
    
    # Clump and harmonize MVMR instrument SNPs across GWAS -------------------
    
    dat <- harmonise_snps(snplist = snplist,
                          t2d_instrument = t2d_instrument,
                          t2d = t2d,
                          feature = feature,
                          outcome = outcome)
    
    # Perform MVMR -----------------------------------------------------------
    
    mvmr_in <- MVMR::format_mvmr(BXGs = dat[,c("beta.t2d","beta.feature")],
                                 BYG = dat[,c("beta.outcome")],
                                 seBXGs = dat[,c("se.t2d","se.feature")],
                                 seBYG = dat[,c("se.outcome")],
                                 RSID = dat[,c("SNP")])
    
    if (sum(mvmr_in$SNP %in% snplist[snplist$instrument=="feature",]$SNP)>0 & sum(mvmr_in$SNP %in% snplist[snplist$instrument=="t2d",]$SNP)>0) {
      
      mvmr_out <- MVMR::ivw_mvmr(r_input = mvmr_in, 
                                 gencov = 0)
      
      mvmr_stren <- MVMR::strength_mvmr(r_input = mvmr_in, 
                                        gencov = 0)
      
      mvmr_plei <- MVMR::pleiotropy_mvmr(r_input = mvmr_in,
                                         gencov = 0)
      
      results[nrow(results)+1,] <- c(analysis,"t2d",j,"direct",mvmr_out[1,c(1,2,4)],paste0(dat$SNP,collapse = ";"),mvmr_plei$Qstat,mvmr_plei$Qpval,mvmr_stren$exposure1)
      results[nrow(results)+1,] <- c(analysis,i,j,"direct",mvmr_out[2,c(1,2,4)],paste0(dat$SNP,collapse = ";"),mvmr_plei$Qstat,mvmr_plei$Qpval,mvmr_stren$exposure2)
      
      # Perform UVMR for t2d on outcome ----------------------------------------
      
      uvmr_in <- dat[dat$SNP %in% snplist[snplist$instrument=="t2d",]$SNP,]
      
      uvmr_in <- MendelianRandomization::mr_input(bx = uvmr_in$beta.t2d,
                                                  bxse = uvmr_in$se.t2d,
                                                  by = uvmr_in$beta.outcome,
                                                  byse = uvmr_in$se.outcome,
                                                  snps = uvmr_in$SNP,
                                                  exposure = t2d_instrument,
                                                  outcome = j,
                                                  effect_allele = uvmr_in$effect_allele,
                                                  other_allele = uvmr_in$other_allele)
      
      uvmr_out <- MendelianRandomization::mr_ivw(uvmr_in)
      
      results[nrow(results)+1,] <- c(analysis,"t2d",j,"total_restricted",uvmr_out@Estimate,uvmr_out@StdError,uvmr_out@Pvalue,paste0(uvmr_in$snps,collapse = ";"),rep(NA,3))
      
      # Perform UVMR for feature on outcome ------------------------------------
      
      uvmr_in <- dat[dat$SNP %in% snplist[snplist$instrument=="feature",]$SNP,]
      
      uvmr_in <- MendelianRandomization::mr_input(bx = uvmr_in$beta.feature,
                                                  bxse = uvmr_in$se.feature,
                                                  by = uvmr_in$beta.outcome,
                                                  byse = uvmr_in$se.outcome,
                                                  snps = uvmr_in$SNP,
                                                  exposure = i,
                                                  outcome = j,
                                                  effect_allele = uvmr_in$effect_allele,
                                                  other_allele = uvmr_in$other_allele)
      
      uvmr_out <- MendelianRandomization::mr_ivw(uvmr_in)
      
      results[nrow(results)+1,] <- c(analysis,i,j,"total_restricted",uvmr_out@Estimate,uvmr_out@StdError,uvmr_out@Pvalue,paste0(uvmr_in$snps,collapse = ";"),rep(NA,3))
      
      # Perform UVMR for t2d on feature ----------------------------------------
      
      uvmr_in <- dat[dat$SNP %in% snplist[snplist$instrument=="t2d",]$SNP,]
      
      uvmr_in <- MendelianRandomization::mr_input(bx = uvmr_in$beta.t2d,
                                                  bxse = uvmr_in$se.t2d,
                                                  by = uvmr_in$beta.feature,
                                                  byse = uvmr_in$se.feature,
                                                  snps = uvmr_in$SNP,
                                                  exposure = t2d_instrument,
                                                  outcome = j,
                                                  effect_allele = uvmr_in$effect_allele,
                                                  other_allele = uvmr_in$other_allele)
      
      uvmr_out <- MendelianRandomization::mr_ivw(uvmr_in)
      
      results[nrow(results)+1,] <- c(analysis,"t2d",i,"total_restricted",uvmr_out@Estimate,uvmr_out@StdError,uvmr_out@Pvalue,paste0(uvmr_in$snps,collapse = ";"),rep(NA,3))
      
      # Perform UVMR for feature on t2d ------------------------------------
      
      uvmr_in <- dat[dat$SNP %in% snplist[snplist$instrument=="feature",]$SNP,]
      
      uvmr_in <- MendelianRandomization::mr_input(bx = uvmr_in$beta.feature,
                                                  bxse = uvmr_in$se.feature,
                                                  by = uvmr_in$beta.t2d,
                                                  byse = uvmr_in$se.t2d,
                                                  snps = uvmr_in$SNP,
                                                  exposure = i,
                                                  outcome = j,
                                                  effect_allele = uvmr_in$effect_allele,
                                                  other_allele = uvmr_in$other_allele)
      
      uvmr_out <- MendelianRandomization::mr_ivw(uvmr_in)
      
      results[nrow(results)+1,] <- c(analysis,i,"t2d","total_restricted",uvmr_out@Estimate,uvmr_out@StdError,uvmr_out@Pvalue,paste0(uvmr_in$snps,collapse = ";"),rep(NA,3))
      
    }
    
  }
  
}

data.table::fwrite(results,"output/mvmr_results.csv",row.names = FALSE)