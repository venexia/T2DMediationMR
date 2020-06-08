rm(list=ls())
graphics.off()

# Specify paths ----------------------------------------------------------------

source("code/specify_paths.R")

# Load instruments -------------------------------------------------------------

instruments <- data.table::fread("data/instruments_all.txt",
                                 stringsAsFactors = FALSE,
                                 data.table = FALSE)

# Generate chromosome position annotation for instruments ----------------------

instruments$chrpos <- paste(instruments$chr, instruments$pos, sep = ":")

# Create results data frame ----------------------------------------------------

results <- NULL

for (j in c("t2d","cad","pad")) {
  
  # Exclude instrument for outcome ---------------------------------------------
  
  df <- instruments[!grepl(j,instruments$exposure),]
  
  if (j=="cad") {
    
    outcome <- TwoSampleMR::extract_outcome_data(snps = df$SNP,
                                                 outcomes = "ieu-a-7",
                                                 proxies = FALSE)
    
    outcome$outcome <- j
    colnames(outcome) <- gsub(".outcome","",colnames(outcome))
    outcome[,c("id","originalname","outcome.deprecated","mr_keep","data_source")] <- NULL
    
  } else {
  
  # Load outcome GWAS ----------------------------------------------------------
  
  outcome <- data.table::fread(paste0("data/gwas-",j,".txt"),
                               stringsAsFactors = FALSE,
                               data.table = FALSE)
  
  # Generate chromosome position annotation for outcome ------------------------
  
  outcome$SNP <- NULL
  outcome$chrpos <- paste(outcome$chr,outcome$pos,sep = ":")
  
  # Restrict outcome GWAS to required SNPs and add rsIDs -----------------------
  
  outcome <- outcome[outcome$chrpos %in% df$chrpos,]
  outcome <- merge(outcome,df[,c("SNP","chrpos")],by=c("chrpos"))
  outcome <- unique(outcome)
  outcome$exposure <- NULL
  outcome$outcome <- j
  
  }
  
  # Perform analysis for each feature-outcome combination ------------------------
  
  for (i in unique(df$exposure)) {
    
    print(paste0("Feature: ",i))
    
    ## Create exposure dataset ---------------------------------------------------
    
    exp <- TwoSampleMR::format_data(dat = df[df$exposure==i,],
                                    type = "exposure",
                                    phenotype_col = "exposure")
    
    ## Note: instruments already clumped
    
    ## Create outcome dataset ----------------------------------------------------
    
    out <- TwoSampleMR::format_data(dat = outcome,
                                    type = "outcome",
                                    snps = exp$SNP,
                                    phenotype_col = "outcome")
    
    ## Harmonise data ----------------------------------------------------------
    
    dat <- TwoSampleMR::harmonise_data(exposure_dat = exp, 
                                       outcome_dat = out)
    
    ## Perform MR --------------------------------------------------------------
    
    tmp_results <- TwoSampleMR::mr(dat = dat)
    
    ## Record SNPs -------------------------------------------------------------
    
    tmp_results$snps_in <- paste(exp$SNP, collapse = ";")
    tmp_results$snps_mr <- paste(dat$SNP, collapse = ";")
    
    ## Save to results data frame ----------------------------------------------
    
    results <- rbind(results,tmp_results)
    
  }
  
}

data.table::fwrite(results,"output/uvmr_results.csv",row.names = FALSE)
