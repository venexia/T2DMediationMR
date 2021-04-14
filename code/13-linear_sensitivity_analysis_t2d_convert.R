rm(list=ls())
graphics.off()

# Specify paths ----------------------------------------------------------------

source("code/specify_paths.R")
source("code/fn-Isq.R")

# Load source data info --------------------------------------------------------

gwas <- data.table::fread("raw/gwas.csv",
                          stringsAsFactors = FALSE,
                          data.table = FALSE)

gwas <- gwas[gwas$source!="exclude_feature" & !grepl("ukb",gwas$ieugwas) & !(gwas$trait %in% c("t2d","pad","cad")),]
rownames(gwas) <- NULL

# Load outcome data ------------------------------------------------------------

outcome <- data.table::fread("raw/t2d_imputed.txt.gz", data.table = FALSE)

ncase <- 24884
ncont <- 437996
mu <- ncase / (ncase + ncont)

outcome$BETA_convert <- outcome$BETA / (mu * (1 - mu))
outcome$SE_convert <- outcome$SE / (mu * (1 - mu))

# Load instruments -------------------------------------------------------------

instruments <- data.table::fread("data/instruments_all.txt", data.table = FALSE)
instruments <- instruments[instruments$exposure %in% gwas$trait,]

results_linear <- NULL
plei_linear <- NULL

for (exposure in unique(instruments$exposure)) {
  
  # Format exposure --------------------------------------------------------------
  
  exp <- TwoSampleMR::format_data(dat = instruments[instruments$exposure==exposure,],
                                  type = "exposure",
                                  phenotype_col = "exposure")
  
  # Clump data ----------------------------------------------------------------
  
  exp <- suppressMessages(TwoSampleMR::clump_data(exp,
                                                  clump_kb = 10000,
                                                  clump_r2 = 0.001,
                                                  pop = "EUR"))
  
  # Format outcome -----------------------------------------------------------
  
  out <- TwoSampleMR::format_data(outcome,
                                  type = "outcome",
                                  snps = exp$SNP,
                                  beta_col = "BETA_convert",
                                  se_col = "SE_convert",
                                  eaf_col = "A1FREQ",
                                  effect_allele_col = "ALLELE1",
                                  other_allele_col = "ALLELE0",
                                  pval_col = "P_LINREG",
                                  chr_col = "CHR",
                                  pos_col = "BP")
  
  out$outcome <- "t2d_linear_convert"
  
  # Harmonise data -----------------------------------------------------------
  
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

# Save -------------------------------------------------------------------------

data.table::fwrite(results_linear,"output/results_linear_convert.csv")
data.table::fwrite(plei_linear,"output/plei_linear_convert.csv")