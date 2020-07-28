rm(list=ls())
graphics.off()

# Specify paths ----------------------------------------------------------------

source("code/specify_paths.R")

# Load instruments -------------------------------------------------------------

instruments <- data.table::fread("data/instruments_all.txt",
                                 stringsAsFactors = FALSE,
                                 data.table = FALSE)

# Specify outcome list ---------------------------------------------------------

finn <- c("finn-a-I_INFECT_PARASIT",
          "finn-a-II_NEOPLASM",
          "finn-a-III_BLOOD_IMMUN",
          "finn-a-IV_ENDOCRIN_NUTRIT",
          "finn-a-V_MENTAL_BEHAV",
          "finn-a-VI_NERVOUS",
          "finn-a-VII_EYE_ADNEXA",
          "finn-a-VIII_EAR_MASTOID",
          "finn-a-IX_CIRCULATORY",
          "finn-a-X_RESPIRATORY",
          "finn-a-XI_DIGESTIVE",
          "finn-a-XII_SKIN_SUBCUTAN",
          "finn-a-XIII_MUSCULOSKELET",
          "finn-a-XIV_GENITOURINARY",
          "finn-a-XV_PREGNANCY_BIRTH",
          "finn-a-XVI_PERINATAL",
          "finn-a-XVII_MALFORMAT_ABNORMAL",
          "finn-a-XVIII_MISCFINDINGS",
          "finn-a-XIX_INJURY_POISON",
          "finn-a-XX_EXTERNAL_MORB_MORT",
          "finn-a-XXI_HEALTHFACTORS")

# Create results data frame ----------------------------------------------------

results <- NULL

# Perform analysis for each feature-outcome combination ------------------------

for (i in unique(instruments$exposure)) {
  
  print(paste0("Feature: ",i))
  
  ## Create exposure dataset -------------------------------------------------
  
  exp <- TwoSampleMR::format_data(dat = instruments[instruments$exposure==i,],
                                  type = "exposure",
                                  phenotype_col = "exposure")
  
  ## Clump data --------------------------------------------------------------
  
  exp <- TwoSampleMR::clump_data(exp,
                                 clump_kb = 10000,
                                 clump_r2 = 0.001,
                                 pop = "EUR")
  
  ## Create outcome dataset --------------------------------------------------
  
  out <- TwoSampleMR::extract_outcome_data(snps = exp$SNP,
                                           outcomes = finn,
                                           proxies = FALSE)
  
  ## Harmonise data ----------------------------------------------------------
  
  dat <- TwoSampleMR::harmonise_data(exposure_dat = exp, 
                                     outcome_dat = out)
  
  ## Perform MR --------------------------------------------------------------
  
  tmp_results <- TwoSampleMR::mr(dat = dat, method_list = "mr_ivw","mr_wald_ratio")
  
  if (nrow(tmp_results) > 0 ) {
    
    ## Record SNPs -------------------------------------------------------------
    
    tmp_results$nsnp_in <- nrow(exp)
    tmp_results$snps_in <- paste(exp$SNP, collapse = ";")
    tmp_results$snps_mr <- paste(dat$SNP, collapse = ";")
    
    ## Save to results data frame ----------------------------------------------
    
    results <- rbind(results,tmp_results)
  
  }
  
}

data.table::fwrite(results,"output/assess_assumptions.csv",row.names = FALSE)
