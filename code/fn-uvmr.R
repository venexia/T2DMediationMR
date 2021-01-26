uvmr <- function(exposure, outcome) {

  # Load source data info ------------------------------------------------------
  
  gwas <- data.table::fread("raw/gwas.csv",
                            stringsAsFactors = FALSE,
                            data.table = FALSE)
  
  # Load instruments -----------------------------------------------------------
  
  instruments <- data.table::fread("data/instruments_all.txt", data.table = FALSE)
  
  ## Create exposure dataset ---------------------------------------------------
  
  exp <- TwoSampleMR::format_data(dat = instruments[instruments$exposure==exposure,],
                                  type = "exposure",
                                  phenotype_col = "exposure")
  
  ## Clump data ----------------------------------------------------------------
  
  exp <- suppressMessages(TwoSampleMR::clump_data(exp,
                                                  clump_kb = 10000,
                                                  clump_r2 = 0.001,
                                                  pop = "EUR"))
  
  # Extract outcome data -------------------------------------------------------
  
  if (is.na(gwas[gwas$trait==outcome,]$ieugwas)) {
    
    out <- data.table::fread(paste0("data/gwas-",outcome,".txt"), data.table = FALSE)
    
    if (sum(grepl("rs",out$SNP))==0) {
      
      out$SNP <- NULL
      tmp <- exp[,c("SNP","chr.exposure","pos.exposure")]
      colnames(tmp) <- c("SNP","chr","pos")
      out <- merge(out, tmp, by = c("chr","pos"))
      
    }
    
    if (nrow(out)<1) {
      
      out <- NULL
      
    } else {
      
      out <- TwoSampleMR::format_data(dat = out,
                                      type = "outcome",
                                      snps = exp$SNP,
                                      phenotype_col = "exposure")
      
    }
    
  } else {
    
    out <- TwoSampleMR::extract_outcome_data(snps = exp$SNP,
                                             outcomes = gwas[gwas$trait==outcome,]$ieugwas,
                                             proxies = FALSE)
    
  }
  
  if (!is.null(out)) {
    
    ## Format outcome data -----------------------------------------------------
    
    out$outcome <- outcome
    
    ## Harmonise data ----------------------------------------------------------
    
    dat <- suppressMessages(TwoSampleMR::harmonise_data(exposure_dat = exp, 
                                                        outcome_dat = out))
    
    ## Perform MR --------------------------------------------------------------
    
    results <- suppressMessages(TwoSampleMR::mr(dat = dat))
    
    ## Record SNPs -------------------------------------------------------------
    
    results$snps_in <- paste(exp$SNP, collapse = ";")
    results$snps_mr <- paste(dat$SNP, collapse = ";")
    
    ## Calculate Isq -----------------------------------------------------------
    
    dat_isq <- exp[exp$SNP %in% dat[dat$mr_keep==TRUE,]$SNP,]
    
    isq <- Isq(dat_isq$beta.exposure,dat_isq$se.exposure)
    
    results$Isq <- isq
    
    # Perform Egger intercept test ---------------------------------------------
    
    plei <- TwoSampleMR::mr_pleiotropy_test(dat)
    
  } else {
    
    results <- NULL
    plei <- NULL
    
  }
  
  r <- list(results, plei)
  
  return(r)
  
}