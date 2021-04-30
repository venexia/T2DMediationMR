rm(list=ls())
graphics.off()

# Specify paths ----------------------------------------------------------------

source("code/specify_paths.R", echo = TRUE)

# Source functions -------------------------------------------------------------

source("code/fn-extract_snps.R", echo = TRUE)
source("code/fn-harmonise_snps.R", echo = TRUE)

# Load instruments -----------------------------------------------------------

instruments_all <- data.table::fread("data/instruments_all.txt",
                                     data.table = FALSE)

# Load filtered features list --------------------------------------------------

features <- data.table::fread("output/evidence_summary.csv", data.table = FALSE)
features <- features[(features$feature_t2d_evidence+features$t2d_feature_evidence)==1 & (features$feature_pad_evidence+features$feature_cad_evidence)>0,]
features <- features[!is.na(features$trait),]

# Specify results data frame ---------------------------------------------------

results <- data.frame(analysis = character(),
                      exposure = character(),
                      outcome = character(),
                      effect = character(),
                      estimate = numeric(),
                      se = numeric(),
                      pval = numeric(),
                      snps = character(),
                      Qstat = numeric(),
                      Qpval = numeric(),
                      Qdf = numeric(),
                      condF = numeric(),
                      stringsAsFactors = FALSE)

# Load source data info --------------------------------------------------------

gwas <- data.table::fread("raw/gwas.csv",
                          stringsAsFactors = FALSE,
                          data.table = FALSE)

for (k in c("t2d","t2d_linear")) {
  
  for (j in c("pad","cad")) {
    
    # Determine features of interest ---------------------------------------------
    
    tmp <- features[,c("trait",paste0("feature_",j,"_evidence"))]
    colnames(tmp) <- c("feature","evidence")
    tmp <- tmp[tmp$evidence==TRUE,]$feature
    
    # Load instruments -----------------------------------------------------------
    
    instruments <- instruments_all[instruments_all$exposure %in% c(tmp,k),]
    
    for (i in tmp) {
      
      analysis <- paste(k, i, j, sep = "/")
      print(paste0("Analysis: ",analysis))
      
      # List MVMR instrument SNPs ----------------------------------------------
      
      snplist <- unique(instruments[instruments$exposure %in% c(i,k),
                                    c("SNP","chrpos")])
      
      snplist$instrument <- ""
      snplist$instrument <- ifelse(snplist$SNP %in% instruments[instruments$exposure==i,]$SNP & snplist$SNP %in% instruments[instruments$exposure==k,]$SNP,"both",snplist$instrument)
      snplist$instrument <- ifelse(snplist$SNP %in% instruments[instruments$exposure==i,]$SNP & !(snplist$SNP %in% instruments[instruments$exposure==k,]$SNP),"feature",snplist$instrument)
      snplist$instrument <- ifelse(!(snplist$SNP %in% instruments[instruments$exposure==i,]$SNP) & snplist$SNP %in% instruments[instruments$exposure==k,]$SNP,k,snplist$instrument)
      
      # Extract MVMR instrument SNPs from GWAS ---------------------------------
      
      t2d <- extract_snps(gwas = k,
                          ieugwas = NA,
                          snplist = snplist,
                          label = k)
      
      colnames(t2d) <- gsub("_linear","",colnames(t2d))
      
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
                            t2d_instrument = k,
                            t2d = t2d,
                            feature = feature,
                            outcome = outcome)
      
      # Perform MVMR -----------------------------------------------------------
      
      mvmr_in <- MVMR::format_mvmr(BXGs = dat[,c("beta.t2d","beta.feature")],
                                   BYG = dat[,c("beta.outcome")],
                                   seBXGs = dat[,c("se.t2d","se.feature")],
                                   seBYG = dat[,c("se.outcome")],
                                   RSID = dat[,c("SNP")])
      
      if (sum(mvmr_in$SNP %in% snplist[snplist$instrument=="feature",]$SNP)>0 & sum(mvmr_in$SNP %in% snplist[snplist$instrument==k,]$SNP)>0) {
        
        mvmr_out <- MVMR::ivw_mvmr(r_input = mvmr_in, 
                                   gencov = 0)
        
        mvmr_stren <- MVMR::strength_mvmr(r_input = mvmr_in, 
                                          gencov = 0)
        
        mvmr_plei <- MVMR::pleiotropy_mvmr(r_input = mvmr_in,
                                           gencov = 0)
        
        Qdf <- gsub(" DF.*","",gsub(".*on ","",capture.output(mvmr_plei <- MVMR::pleiotropy_mvmr(r_input = mvmr_in,gencov = 0))[2]))
        
        results[nrow(results)+1,] <- c(analysis,k,j,"direct",mvmr_out[1,c(1,2,4)],paste0(dat$SNP,collapse = ";"),mvmr_plei$Qstat,mvmr_plei$Qpval,Qdf,mvmr_stren$exposure1)
        results[nrow(results)+1,] <- c(analysis,i,j,"direct",mvmr_out[2,c(1,2,4)],paste0(dat$SNP,collapse = ";"),mvmr_plei$Qstat,mvmr_plei$Qpval,Qdf,mvmr_stren$exposure2)
        
        # Perform UVMR for t2d on outcome ----------------------------------------
        
        uvmr_in <- dat[dat$SNP %in% snplist[snplist$instrument==k,]$SNP,]
        
        uvmr_in <- MendelianRandomization::mr_input(bx = uvmr_in$beta.t2d,
                                                    bxse = uvmr_in$se.t2d,
                                                    by = uvmr_in$beta.outcome,
                                                    byse = uvmr_in$se.outcome,
                                                    snps = uvmr_in$SNP,
                                                    exposure = "t2d",
                                                    outcome = j,
                                                    effect_allele = uvmr_in$effect_allele,
                                                    other_allele = uvmr_in$other_allele)
        
        uvmr_out <- MendelianRandomization::mr_ivw(uvmr_in)
        
        results[nrow(results)+1,] <- c(analysis,k,j,"total_restricted",uvmr_out@Estimate,uvmr_out@StdError,uvmr_out@Pvalue,paste0(uvmr_in$snps,collapse = ";"),rep(NA,4))
        
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
        
        results[nrow(results)+1,] <- c(analysis,i,j,"total_restricted",uvmr_out@Estimate,uvmr_out@StdError,uvmr_out@Pvalue,paste0(uvmr_in$snps,collapse = ";"),rep(NA,4))
        
        # Perform UVMR for t2d on feature ----------------------------------------
        
        uvmr_in <- dat[dat$SNP %in% snplist[snplist$instrument==k,]$SNP,]
        
        uvmr_in <- MendelianRandomization::mr_input(bx = uvmr_in$beta.t2d,
                                                    bxse = uvmr_in$se.t2d,
                                                    by = uvmr_in$beta.feature,
                                                    byse = uvmr_in$se.feature,
                                                    snps = uvmr_in$SNP,
                                                    exposure = k,
                                                    outcome = j,
                                                    effect_allele = uvmr_in$effect_allele,
                                                    other_allele = uvmr_in$other_allele)
        
        uvmr_out <- MendelianRandomization::mr_ivw(uvmr_in)
        
        results[nrow(results)+1,] <- c(analysis,k,i,"total_restricted",uvmr_out@Estimate,uvmr_out@StdError,uvmr_out@Pvalue,paste0(uvmr_in$snps,collapse = ";"),rep(NA,4))
        
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
        
        results[nrow(results)+1,] <- c(analysis,i,k,"total_restricted",uvmr_out@Estimate,uvmr_out@StdError,uvmr_out@Pvalue,paste0(uvmr_in$snps,collapse = ";"),rep(NA,4))
        
      }
      
    }
    
  }
  
}

data.table::fwrite(results,"output/mvmr_results.csv",row.names = FALSE)