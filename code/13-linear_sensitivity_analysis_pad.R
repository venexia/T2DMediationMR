rm(list=ls())
graphics.off()

# Specify paths ----------------------------------------------------------------

source("code/specify_paths.R")
source("code/fn-Isq.R")

# Load source data info --------------------------------------------------------

outcome <- data.table::fread("raw/PAD.glm.linear", data.table = FALSE)
outcome$chrpos <- paste0(outcome$CHR,":",outcome$POS)
outcome$A0 <- NA
outcome$A0 <- ifelse(outcome$A1==outcome$REF,outcome$ALT,outcome$A0)
outcome$A0 <- ifelse(outcome$A1==outcome$ALT,outcome$REF,outcome$A0)

# Load instruments -------------------------------------------------------------

instruments <- data.table::fread("data/instruments_all.txt", data.table = FALSE)

results_linear <- NULL
plei_linear <- NULL

for (exposure in unique(instruments$exposure)) {
  
  if (!(exposure %in% c("pad","cad"))) {
    
    # Format exposure --------------------------------------------------------------
    
    exp <- TwoSampleMR::format_data(dat = instruments[instruments$exposure==exposure,],
                                    type = "exposure",
                                    phenotype_col = "exposure")
    
    # Clump data ----------------------------------------------------------------
    
    exp <- suppressMessages(TwoSampleMR::clump_data(exp,
                                                    clump_kb = 10000,
                                                    clump_r2 = 0.001,
                                                    pop = "EUR"))
    
    # Format outcome
    
    tmp <- instruments[instruments$SNP %in% exp$SNP & instruments$exposure==exposure,c("SNP","chrpos")]
    out <- merge(outcome, tmp, by = "chrpos")
    
    out <- TwoSampleMR::format_data(out,
                                    type = "outcome",
                                    snps = exp$SNP,
                                    beta_col = "BETA",
                                    se_col = "SE",
                                    eaf_col = "A1_FREQ",
                                    effect_allele_col = "A1",
                                    other_allele_col = "A0",
                                    pval_col = "P",
                                    chr_col = "CHR",
                                    pos_col = "POS")
    
    out$outcome <- "pad_linear"
    
    # Harmonise data ----------------------------------------------------------
    
    dat <- suppressMessages(TwoSampleMR::harmonise_data(exposure_dat = exp, 
                                                        outcome_dat = out))
    
    # Perform MR --------------------------------------------------------------
    
    results <- suppressMessages(TwoSampleMR::mr(dat = dat))
    results$nsnp_ins <- nrow(exp)
    
    # Record SNPs --------------------------------------------------------------
    
    results$snps_in <- paste(exp$SNP, collapse = ";")
    results$snps_mr <- paste(dat$SNP, collapse = ";")
    
    # Calculate Isq ------------------------------------------------------------
    
    dat_isq <- exp[exp$SNP %in% dat[dat$mr_keep==TRUE,]$SNP,]
    
    isq <- Isq(dat_isq$beta.exposure,dat_isq$se.exposure)
    
    results$Isq <- isq
    
    results_linear <- rbind(results, results_linear)
    
    # Perform Egger intercept test ---------------------------------------------
    
    plei <- TwoSampleMR::mr_pleiotropy_test(dat)
    
    plei_linear <- rbind(plei, plei_linear)
    
  }
}

# Save -------------------------------------------------------------------------

data.table::fwrite(results_linear,"output/results_linear.csv")
data.table::fwrite(plei_linear,"output/plei_linear.csv")