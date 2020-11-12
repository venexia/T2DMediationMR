rm(list=ls())
graphics.off()

# Specify paths ----------------------------------------------------------------

source("code/specify_paths.R")
source("code/fn-Isq.R")

# Load instruments -------------------------------------------------------------

instruments <- data.table::fread("data/instruments_all.txt", data.table = FALSE)

# # Extract instruments from outcome data ----------------------------------------
# 
# cad <- TwoSampleMR::extract_outcome_data(snps = instruments$SNP,
#                                          outcomes = "ieu-a-7",
#                                          proxies = FALSE)
# 
# cad$outcome <- "cad"
# colnames(cad) <- gsub(".outcome","",colnames(cad))
# cad[,c("id","originalname","outcome.deprecated","mr_keep","data_source")] <- NULL
# 
# data.table::fwrite(cad, "data/outcome_cad.txt")
#   
# for (x in c("t2d","pad")) {
#   
#   # Load outcome GWAS ----------------------------------------------------------
#   
#   outcome <- data.table::fread(paste0("data/gwas-",x,".txt"), data.table = FALSE)
#   
#   # Generate chromosome position annotation for outcome ------------------------
#   
#   outcome$SNP <- NULL
#   outcome$chrpos <- paste(outcome$chr,outcome$pos,sep = ":")
#   
#   # Restrict outcome GWAS to required SNPs and add rsIDs -----------------------
#   
#   outcome <- outcome[outcome$chrpos %in% instruments$chrpos,]
#   outcome <- merge(outcome,instruments[,c("SNP","chrpos")],by=c("chrpos"))
#   outcome <- unique(outcome)
#   outcome$exposure <- NULL
#   outcome$outcome <- x
# 
#   # Name data ------------------------------------------------------------------
#  
#   data.table::fwrite(outcome, paste0("data/outcome_",x,".txt")) 
#   assign(x, outcome)
# 
# }

t2d <- data.table::fread("data/outcome_t2d.txt", data.table = FALSE)
pad <- data.table::fread("data/outcome_pad.txt", data.table = FALSE)
cad <- data.table::fread("data/outcome_cad.txt", data.table = FALSE)

# Create results data frame ----------------------------------------------------

run <- data.frame(exposure = rep(unique(instruments$exposure), each = 3),
                  outcome = rep(c("pad","t2d","cad"), times = length(unique(instruments$exposure))),
                  stringsAsFactors = FALSE)

# results <- NULL
# plei <- NULL

results <- data.table::fread("output/uvmr_results.csv", data.table = FALSE)
plei <- data.table::fread("output/uvmr_pleiotropy.csv", data.table = FALSE)

# 86

for (k in 86) {
  
  print(paste0(k," of ",nrow(run)))

  i <- run$exposure[k]
  j <- run$outcome[k]
  outcome <- get(j)
  
  ## Create exposure dataset -------------------------------------------------

  exp <- TwoSampleMR::format_data(dat = instruments[instruments$exposure==i,],
                                  type = "exposure",
                                  phenotype_col = "exposure")
  
  ## Clump data --------------------------------------------------------------
  
  exp <- suppressMessages(TwoSampleMR::clump_data(exp,
                                 clump_kb = 10000,
                                 clump_r2 = 0.001,
                                 pop = "EUR"))
  # Extract outcome data -------------------------------------------------------
  
  out <- TwoSampleMR::format_data(dat = outcome,
                                  type = "outcome",
                                  snps = exp$SNP,
                                  phenotype_col = "outcome")
  
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
  
  ## Save to results data frame ----------------------------------------------
  
  results <- rbind(results, tmp_results)
  plei <- rbind(plei, tmp_plei)
  
}

data.table::fwrite(results,"output/uvmr_results.csv",row.names = FALSE)
data.table::fwrite(plei,"output/uvmr_pleiotropy.csv",row.names = FALSE)
