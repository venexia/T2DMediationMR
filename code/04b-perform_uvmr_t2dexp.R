rm(list=ls())
graphics.off()

# Specify paths ----------------------------------------------------------------

source("code/specify_paths.R")
source("code/fn-Isq.R")

# Load source data info --------------------------------------------------------

gwas <- data.table::fread("raw/gwas.csv",
                          stringsAsFactors = FALSE,
                          data.table = FALSE)

gwas <- gwas[!(gwas$trait %in% c("t2d","pad","cad","t2d_udler")),]

# Load instruments -------------------------------------------------------------

instruments <- data.table::fread("data/instruments_all.txt", data.table = FALSE)

# Create results data frame ----------------------------------------------------

results <- data.table::fread("output/uvmr_results.csv", data.table = FALSE)
plei <- data.table::fread("output/uvmr_pleiotropy.csv", data.table = FALSE)

## Create exposure dataset -------------------------------------------------

exp <- TwoSampleMR::format_data(dat = instruments[instruments$exposure=="t2d",],
                                type = "exposure",
                                phenotype_col = "exposure")

## Clump data --------------------------------------------------------------

exp <- suppressMessages(TwoSampleMR::clump_data(exp,
                                                clump_kb = 10000,
                                                clump_r2 = 0.001,
                                                pop = "EUR"))

# Create snplist ---------------------------------------------------------------

snplist <- instruments[instruments$exposure=="t2d" & instruments$SNP %in% exp$SNP,c("SNP","chrpos")]

# Run for each outcome ---------------------------------------------------------

for (k in 1:nrow(gwas)) {
  
  print(paste0(k," of ",nrow(gwas)))
  
  # Extract outcome data -------------------------------------------------------
  
  ieugwas <- gwas$ieugwas[k]
  trait <- gwas$trait[k]
    
  if (ieugwas=="") {
    
    out <- data.table::fread(paste0("data/gwas-",trait,".txt"),
                            stringsAsFactors = FALSE,
                            data.table = FALSE)
    
    if (sum(grepl("rs",out$SNP))==0) {
      
      out$chrpos <- paste(out$chr, out$pos, sep = ":")
      
      out <- out[out$chrpos %in% snplist$chrpos,]
      
      out$SNP <- NULL
      
      out <- merge(out, snplist, by = "chrpos")
      
    } else {
      
      out <- out[out$SNP %in% snplist$SNP,]
      
    }
    
    out <- TwoSampleMR::format_data(dat = out,
                                    type = "outcome",
                                    snps = exp$SNP,
                                    phenotype_col = "outcome")
    
  } else {
    
    out <- TwoSampleMR::extract_outcome_data(outcomes = ieugwas,
                                              snps = exp$SNP,
                                              proxies = FALSE)
  }
 
  if (!is.null(out)) {
    
  out$outcome <- trait
  
  ## Harmonise data ------------------------------------------------------------
  
  dat <- suppressMessages(TwoSampleMR::harmonise_data(exposure_dat = exp, 
                                                      outcome_dat = out))
  
  ## Perform MR --------------------------------------------------------------
  
  tmp_results <- suppressMessages(TwoSampleMR::mr(dat = dat))
  
  ## Record SNPs -------------------------------------------------------------
  
  tmp_results$snps_in <- paste(exp$SNP, collapse = ";")
  tmp_results$snps_mr <- paste(dat$SNP, collapse = ";")
  
  ## Calculate Isq -----------------------------------------------------------
  
  dat_isq <- exp[exp$SNP %in% dat[dat$mr_keep==TRUE,]$SNP,]
  
  isq <- Isq(dat_isq$beta.exposure,dat_isq$se.exposure)
  
  tmp_results$Isq <- isq
  
  # Perform Egger intercept test ---------------------------------------------
  
  tmp_plei <- TwoSampleMR::mr_pleiotropy_test(dat)
  
  } else {
    
    tmp_results <- data.frame(id.exposure = "", 
                             id.outcome = "", 
                             outcome = trait,
                             exposure = "t2d",
                             method = "",
                             nsnp = NA,
                             b = NA,
                             se = NA,
                             pval = NA,
                             snps_in = "",
                             snps_mr = "",
                             Isq = NA,
                             stringsAsFactors = FALSE)
    
    tmp_plei <- data.frame(id.exposure = "", 
                           id.outcome = "", 
                           outcome = trait,
                           exposure = "t2d", 
                           egger_intercept = NA, 
                           se = NA, 
                           pval = NA, 
                           stringsAsFactors = FALSE)
    
  }
  
  ## Save to results data frame ----------------------------------------------
  
  results <- rbind(results, tmp_results)
  plei <- rbind(plei, tmp_plei)
  
}

data.table::fwrite(results,"output/uvmr_results.csv",row.names = FALSE)
data.table::fwrite(plei,"output/uvmr_pleiotropy.csv",row.names = FALSE)
