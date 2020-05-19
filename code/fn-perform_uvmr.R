perform_uvmr <- function(analysis,exposure,outcome,outcome_label) {
  
  # Determine whether the exposure and outcome is a feature
  
  if (exposure %in% c("t2d","outcome")) {
    x <- exposure
  } else {
    x <- "feature"
  }
  
  if (outcome %in% c("t2d","outcome")) {
    y <- outcome
  } else {
    y <- "feature"
  }
  
  # Obtain instrument for exposure
  
  instrument <- get(paste0("instrument_",x))
  
  # Format data for MR
  
  uvmr_in <- MendelianRandomization::mr_input(bx = df[df$SNP %in% instrument,c(paste0("beta.",x))],
                                              bxse = df[df$SNP %in% instrument,c(paste0("se.",x))],
                                              by = df[df$SNP %in% instrument,c(paste0("beta.",y))],
                                              byse = df[df$SNP %in% instrument,c(paste0("se.",y))])
  
  # Perform MR
  
  uvmr_out <- MendelianRandomization::mr_ivw(uvmr_in)
  
  # Return results
  
  results <- c(analysis,exposure,outcome_label,"total",uvmr_out@Estimate,uvmr_out@StdError,uvmr_out@Pvalue,paste0(instrument,collapse = ";"),rep(NA,2))
  
  return(results)
  
}