mvmr <- function(exposure, mediator, outcome) {
  
  # Define analysis ------------------------------------------------------------
  
  analysis <- paste(exposure, mediator, outcome, sep = "/")
  
  # Extract exposure instrument ------------------------------------------------
  
  exposure_ins <- NULL
  
  if (exposure %in% c("t2d","t2d_linear")) {
    
    exposure_ins <- data.table::fread(paste0("data/gwas-",exposure,".csv.gz"), data.table = FALSE)
    exposure_ins <- exposure_ins[exposure_ins$pval<5e-8,]
    
  } else {
    
    exposure_ins <- TwoSampleMR::extract_instruments(outcomes = exposure,
                                                     clump = FALSE)
    colnames(exposure_ins) <- gsub(".exposure","",colnames(exposure_ins))
    
  }
  
  # Extract mediator instrument ------------------------------------------------
  
  mediator_ins <- NULL
  
  if (mediator %in% c("t2d","t2d_linear")) {
    
    mediator_ins <- data.table::fread(paste0("data/gwas-",mediator,".csv.gz"), data.table = FALSE)
    mediator_ins <- mediator_ins[mediator_ins$pval<5e-8,]
    
  } else {
    
    mediator_ins <- TwoSampleMR::extract_instruments(outcomes = mediator,
                                                     clump = FALSE)
    colnames(mediator_ins) <- gsub(".exposure","",colnames(mediator_ins))
  }
  
  # Define snplist for MVMR ----------------------------------------------------
  
  snplist <- data.frame(SNP = unique(c(exposure_ins$SNP,mediator_ins$SNP)), 
                        instrument = "", 
                        stringsAsFactors = FALSE)
  
  snplist$instrument <- ifelse(snplist$SNP %in% exposure_ins$SNP & snplist$SNP %in% mediator_ins$SNP,
                               "both", snplist$instrument)
  
  snplist$instrument <- ifelse(snplist$SNP %in% exposure_ins$SNP & !(snplist$SNP %in% mediator_ins$SNP),
                               "exposure", snplist$instrument)
  
  snplist$instrument <- ifelse(!(snplist$SNP %in% exposure_ins$SNP) & snplist$SNP %in% mediator_ins$SNP,
                               "mediator", snplist$instrument)
  
  # Extract exposure SNP list --------------------------------------------------
  
  exposure_gwas <- NULL
  
  if (exposure %in% c("t2d","t2d_linear")) {
    
    exposure_gwas <- data.table::fread(paste0("data/gwas-",exposure,".csv.gz"), data.table = FALSE)
    exposure_gwas <- exposure_gwas[exposure_gwas$SNP %in% snplist$SNP,]
    
  } else {
    
    exposure_gwas <- TwoSampleMR::extract_outcome_data(snps = snplist$SNP,
                                                       outcomes = exposure,
                                                       proxies = TRUE)
    
    colnames(exposure_gwas) <- gsub(".outcome","",colnames(exposure_gwas))
    
  }
  
  # Extract mediator SNP list --------------------------------------------------
  
  mediator_gwas <- NULL
  
  if (mediator %in% c("t2d","t2d_linear")) {
    
    mediator_gwas <- data.table::fread(paste0("data/gwas-",mediator,".csv.gz"), data.table = FALSE)
    mediator_gwas <- mediator_gwas[mediator_gwas$SNP %in% snplist$SNP,]
    
  } else {
    
    mediator_gwas <- TwoSampleMR::extract_outcome_data(snps = snplist$SNP,
                                                       outcomes = mediator,
                                                       proxies = TRUE)
    
    colnames(mediator_gwas) <- gsub(".outcome","",colnames(mediator_gwas))
    
  }
  
  # Clump SNPs against p-values from smaller instrument ------------------------
  
  clumped <- NULL
  
  if (nrow(exposure_ins)<=nrow(mediator_ins)) {
    clumped <- exposure_gwas[,c("SNP","pval")]
    colnames(clumped) <- c("rsid","pval")
  } 
  
  if (nrow(exposure_ins)>nrow(mediator_ins)) {
    clumped <- mediator_gwas[,c("SNP","pval")]
    colnames(clumped) <- c("rsid","pval")
  }
  
  clumped <- na.omit(clumped)
  
  clumped <- ieugwasr::ld_clump(clumped,
                                clump_kb = 10000,
                                clump_r2 = 0.001)
  
  # Restrict exposure and mediator GWAS to clumped SNPs ------------------------
  
  exposure_gwas <- exposure_gwas[exposure_gwas$SNP %in% clumped$rsid,]
  mediator_gwas <- mediator_gwas[mediator_gwas$SNP %in% clumped$rsid,]
  
  # Extract outcome SNP list ---------------------------------------------------
  
  outcome_gwas <- data.table::fread(paste0("data/gwas-",outcome,".csv.gz"), data.table = FALSE)
  outcome_gwas <- unique(outcome_gwas[outcome_gwas$SNP %in% clumped$rsid,])
  
  # Make analysis dataset ------------------------------------------------------
  
  dat <- exposure_gwas[,c("SNP","effect_allele","other_allele","beta","se","pval")]
  colnames(dat) <- c("SNP","effect_allele.exposure","other_allele.exposure","beta.exposure","se.exposure","pval.exposure")
  
  dat <- merge(dat, mediator_gwas[,c("SNP","effect_allele","other_allele","beta","se","pval")], by = "SNP")
  colnames(dat) <- c("SNP",
                     "effect_allele.exposure","other_allele.exposure","beta.exposure","se.exposure","pval.exposure",
                     "orig.effect_allele.mediator","orig.other_allele.mediator","orig.beta.mediator","se.mediator","pval.mediator")
  
  dat <- merge(dat, outcome_gwas[,c("SNP","effect_allele","other_allele","beta","se","pval")], by = "SNP")
  colnames(dat) <- c("SNP",
                     "effect_allele","other_allele","beta.exposure","se.exposure","pval.exposure",
                     "orig.effect_allele.mediator","orig.other_allele.mediator","orig.beta.mediator","se.mediator","pval.mediator",
                     "orig.effect_allele.outcome","orig.other_allele.outcome","orig.beta.outcome","se.outcome","pval.outcome")
  
  # Make all effects exposure increasing --------------------------------------------
  
  dat$beta.mediator <- ifelse(dat$effect_allele==dat$orig.effect_allele.mediator & dat$other_allele==dat$orig.other_allele.mediator,
                              dat$orig.beta.mediator,
                              -1*dat$orig.beta.mediator)
  
  dat$beta.outcome <- ifelse(dat$effect_allele==dat$orig.effect_allele.outcome & dat$other_allele==dat$orig.other_allele.outcome,
                             dat$orig.beta.outcome,
                             -1*dat$orig.beta.outcome)
  
  dat <- dat[,c("SNP","effect_allele","other_allele","beta.exposure","se.exposure","beta.mediator","se.mediator","beta.outcome","se.outcome")]
  
  # Create empty results table -------------------------------------------------
  
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
  
  # Perform MVMR ---------------------------------------------------------------
  
  mvmr_in <- MVMR::format_mvmr(BXGs = dat[,c("beta.exposure","beta.mediator")],
                               BYG = dat[,c("beta.outcome")],
                               seBXGs = dat[,c("se.exposure","se.mediator")],
                               seBYG = dat[,c("se.outcome")],
                               RSID = dat[,c("SNP")])
  
  tmp <- merge(mvmr_in, snplist, by = "SNP")
  
  if (length(unique(tmp$instrument))==3) {
    
    mvmr_out <- MVMR::ivw_mvmr(r_input = mvmr_in, 
                               gencov = 0)
    
    mvmr_stren <- MVMR::strength_mvmr(r_input = mvmr_in, 
                                      gencov = 0)
    
    mvmr_plei <- MVMR::pleiotropy_mvmr(r_input = mvmr_in,
                                       gencov = 0)
    
    Qdf <- gsub(" DF.*","",gsub(".*on ","",capture.output(mvmr_plei <- MVMR::pleiotropy_mvmr(r_input = mvmr_in,gencov = 0))[2]))
    
    results[nrow(results)+1,] <- c(analysis,exposure,outcome,"direct",mvmr_out[1,c(1,2,4)],paste0(dat$SNP,collapse = ";"),mvmr_plei$Qstat,mvmr_plei$Qpval,Qdf,mvmr_stren$exposure1)
    results[nrow(results)+1,] <- c(analysis,mediator,outcome,"direct",mvmr_out[2,c(1,2,4)],paste0(dat$SNP,collapse = ";"),mvmr_plei$Qstat,mvmr_plei$Qpval,Qdf,mvmr_stren$exposure2)
    
    # Perform UVMR for exposure on outcome ----------------------------------------
    
    uvmr_in <- dat[dat$SNP %in% snplist[snplist$instrument!="mediator",]$SNP,]
    
    uvmr_in <- MendelianRandomization::mr_input(bx = uvmr_in$beta.exposure,
                                                bxse = uvmr_in$se.exposure,
                                                by = uvmr_in$beta.outcome,
                                                byse = uvmr_in$se.outcome,
                                                snps = uvmr_in$SNP,
                                                exposure = exposure,
                                                outcome = outcome,
                                                effect_allele = uvmr_in$effect_allele,
                                                other_allele = uvmr_in$other_allele)
    
    uvmr_out <- MendelianRandomization::mr_ivw(uvmr_in)
    
    results[nrow(results)+1,] <- c(analysis,exposure,outcome,"total_restricted",uvmr_out@Estimate,uvmr_out@StdError,uvmr_out@Pvalue,paste0(uvmr_in$snps,collapse = ";"),rep(NA,4))
    
    # Perform UVMR for mediator on outcome ------------------------------------
    
    uvmr_in <- dat[dat$SNP %in% snplist[snplist$instrument!="exposure",]$SNP,]
    
    uvmr_in <- MendelianRandomization::mr_input(bx = uvmr_in$beta.mediator,
                                                bxse = uvmr_in$se.mediator,
                                                by = uvmr_in$beta.outcome,
                                                byse = uvmr_in$se.outcome,
                                                snps = uvmr_in$SNP,
                                                exposure = mediator,
                                                outcome = outcome,
                                                effect_allele = uvmr_in$effect_allele,
                                                other_allele = uvmr_in$other_allele)
    
    uvmr_out <- MendelianRandomization::mr_ivw(uvmr_in)
    
    results[nrow(results)+1,] <- c(analysis,mediator,outcome,"total_restricted",uvmr_out@Estimate,uvmr_out@StdError,uvmr_out@Pvalue,paste0(uvmr_in$snps,collapse = ";"),rep(NA,4))
    
    # Perform UVMR for exposure on mediator ------------------------------------
    
    uvmr_in <- dat[dat$SNP %in% snplist[snplist$instrument!="mediator",]$SNP,]
    
    uvmr_in <- MendelianRandomization::mr_input(bx = uvmr_in$beta.exposure,
                                                bxse = uvmr_in$se.exposure,
                                                by = uvmr_in$beta.mediator,
                                                byse = uvmr_in$se.mediator,
                                                snps = uvmr_in$SNP,
                                                exposure = exposure,
                                                outcome = mediator,
                                                effect_allele = uvmr_in$effect_allele,
                                                other_allele = uvmr_in$other_allele)
    
    uvmr_out <- MendelianRandomization::mr_ivw(uvmr_in)
    
    results[nrow(results)+1,] <- c(analysis,exposure,mediator,"total_restricted",uvmr_out@Estimate,uvmr_out@StdError,uvmr_out@Pvalue,paste0(uvmr_in$snps,collapse = ";"),rep(NA,4))
    
    # Perform UVMR for mediator on exposure ------------------------------------
    
    uvmr_in <- dat[dat$SNP %in% snplist[snplist$instrument!=exposure,]$SNP,]
    
    uvmr_in <- MendelianRandomization::mr_input(bx = uvmr_in$beta.mediator,
                                                bxse = uvmr_in$se.mediator,
                                                by = uvmr_in$beta.exposure,
                                                byse = uvmr_in$se.exposure,
                                                snps = uvmr_in$SNP,
                                                exposure = mediator,
                                                outcome = exposure,
                                                effect_allele = uvmr_in$effect_allele,
                                                other_allele = uvmr_in$other_allele)
    
    uvmr_out <- MendelianRandomization::mr_ivw(uvmr_in)
    
    results[nrow(results)+1,] <- c(analysis,mediator,exposure,"total_restricted",uvmr_out@Estimate,uvmr_out@StdError,uvmr_out@Pvalue,paste0(uvmr_in$snps,collapse = ";"),rep(NA,4))
    
  } 
  
  if (nrow(results)==0) {

    results[nrow(results)+1,] <- c(analysis,exposure,outcome,"direct",rep(NA,3),"",rep(NA,4))
    results[nrow(results)+1,] <- c(analysis,mediator,outcome,"direct",rep(NA,3),"",rep(NA,4))
    results[nrow(results)+1,] <- c(analysis,exposure,outcome,"total_restricted",rep(NA,3),"",rep(NA,4))
    results[nrow(results)+1,] <- c(analysis,mediator,outcome,"total_restricted",rep(NA,3),"",rep(NA,4))
    results[nrow(results)+1,] <- c(analysis,exposure,mediator,"total_restricted",rep(NA,3),"",rep(NA,4))
    results[nrow(results)+1,] <- c(analysis,mediator,exposure,"total_restricted",rep(NA,3),"",rep(NA,4))
    
  }
  
  # Return results -------------------------------------------------------------
  
  return(results)
  
}