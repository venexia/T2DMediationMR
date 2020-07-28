rm(list=ls())
graphics.off()

# Set p-value threshold for features to progress to next analysis --------------

threshold <- 0.05

# Load univariate results ------------------------------------------------------

uvmr <- data.table::fread("output/uvmr_results.csv", data.table = FALSE)

# Load instruments -------------------------------------------------------------

instruments <- data.table::fread("data/instruments_all.txt", data.table = FALSE)

for (i in c("cad","pad")) {

  # Restrict to relavant methods, exposures and outcomes -----------------------

  df <- uvmr[uvmr$method %in% c("Inverse variance weighted","Wald ratio") &
           !(uvmr$exposure %in% c("t2d","t2d_udler","cad","pad")) &
           uvmr$outcome %in% c("t2d",i),
         c("exposure","outcome","pval")]

  # Make outcome name generic --------------------------------------------------
  
  df$outcome <- ifelse(df$outcome==i,"outcome",df$outcome)
  
  # Convert data to wide -------------------------------------------------------
  
  df <- tidyr::pivot_wider(df, 
                           names_from = "outcome", 
                           values_from = c("pval"))


  # Apply thresholds -----------------------------------------------------------
  
  df <- df[df$t2d<threshold & df$outcome<threshold,]
  
  # Restrict instruments -------------------------------------------------------
  
  instruments_restricted <- instruments[instruments$exposure %in% c(df$exposure,"t2d","t2d_udler","cad","pad"),]

  # Save -----------------------------------------------------------------------
  
  data.table::fwrite(instruments_restricted,paste0("data/instruments_",i,".txt"))

}