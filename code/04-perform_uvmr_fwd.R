rm(list=ls())
graphics.off()

# Specify paths ----------------------------------------------------------------

source("code/specify_paths.R")
source("code/fn-Isq.R")
source("code/fn-uvmr.R")

# Load instruments -------------------------------------------------------------

exp <- data.table::fread("data/instruments.csv", data.table = FALSE)

exp <- TwoSampleMR::format_data(dat = exp,
                                type = "exposure",
                                phenotype_col = "exposure",
                                snp_col = "SNP",
                                beta_col = "beta",
                                se_col = "se",
                                eaf_col = "eaf",
                                effect_allele_col = "effect_allele",
                                other_allele_col = "other_allele",
                                pval_col = "pval",
                                units_col = "units",
                                ncase_col = "ncase",
                                ncontrol_col = "ncontrol",
                                samplesize_col = "samplesize",
                                chr_col = "chr",
                                pos_col = "pos")

# Create empty results datasets ------------------------------------------------

results_fwd <- NULL #data.table::fread("output/results.csv", data.table = FALSE)
plei_fwd <- NULL #data.table::fread("output/plei.csv", data.table = FALSE)

for (outcome in c("pad","cad","t2d_linear","t2d")) {
  
  # Load outcome data ----------------------------------------------------------
  
  out <- data.table::fread(paste0("data/gwas-",outcome,".csv.gz"), data.table = FALSE)
  
  # Annotate outcome SNPs with rsIDs -------------------------------------------
  
  if (sum(grepl("rs",out$SNP))==0) {
    out$SNP <- NULL
    tmp <- exp[,c("SNP","chr.exposure","pos.exposure")]
    colnames(tmp) <- c("SNP","chr","pos")
    out <- merge(out, tmp, by = c("chr","pos"))
  }
  
  # Format outcome data --------------------------------------------------------
  
  out <- TwoSampleMR::format_data(dat = out,
                                  type = "outcome",
                                  snps = exp$SNP,
                                  phenotype_col = "exposure")
  
  # Harmonise data -------------------------------------------------------------
  
  dat <- TwoSampleMR::harmonise_data(exposure_dat = exp, outcome_dat = out)
  dat <- dat[dat$exposure!=dat$outcome,]
  
  # Perform MR -----------------------------------------------------------------
  
  results <- TwoSampleMR::mr(dat = dat)
  
  # Calculate Isq --------------------------------------------------------------
  
  dat_isq <- exp[exp$SNP %in% dat[dat$mr_keep==TRUE,]$SNP,]
  
  isq <- Isq(dat_isq$beta.exposure,dat_isq$se.exposure)
  
  results$Isq <- isq
  
  # Perform Egger intercept test -----------------------------------------------
  
  plei <- TwoSampleMR::mr_pleiotropy_test(dat)
  
  # Record results -------------------------------------------------------------
  
  results_fwd <- rbind(results_fwd, results)
  plei_fwd <- rbind(plei_fwd, plei)
  
}

# Save -------------------------------------------------------------------------

data.table::fwrite(results_fwd, "output/results_fwd.csv")
data.table::fwrite(plei_fwd, "output/plei_fwd.csv")