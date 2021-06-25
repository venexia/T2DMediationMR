rm(list=ls())
graphics.off()

# Specify paths ----------------------------------------------------------------

source("code/specify_paths.R")
source("code/fn-Isq.R")
source("code/fn-uvmr.R")

# Load risk factor information -------------------------------------------------

risk_factors <- data.table::fread("data/risk_factors.csv", data.table = FALSE)

# Load exposure data -----------------------------------------------------------

exp <- data.table::fread("data/instruments.csv", data.table = FALSE)

exp <- exp[exp$exposure %in% c("t2d","t2d_linear"),]

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

results_bkwd <- NULL
plei_bkwd <- NULL

# Extract outcome data ---------------------------------------------------------

out <- TwoSampleMR::extract_outcome_data(snps = exp$SNP,
                                         outcomes = risk_factors[risk_factors$id!=risk_factors$trait,]$id,
                                         proxies = FALSE)

out$outcome <- out$id.outcome

# Harmonise data ---------------------------------------------------------------

dat <- TwoSampleMR::harmonise_data(exposure_dat = exp, outcome_dat = out)

# Perform MR -------------------------------------------------------------------

results <- TwoSampleMR::mr(dat = dat)

# Calculate Isq ----------------------------------------------------------------

dat_isq <- exp[exp$SNP %in% dat[dat$mr_keep==TRUE,]$SNP,]
isq <- Isq(dat_isq$beta.exposure,dat_isq$se.exposure)
results$Isq <- isq

# Perform Egger intercept test -------------------------------------------------

plei <- TwoSampleMR::mr_pleiotropy_test(dat)

# Save -------------------------------------------------------------------------

data.table::fwrite(results, "output/results_bkwd.csv")
data.table::fwrite(plei, "output/plei_bkwd.csv")