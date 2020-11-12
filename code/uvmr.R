rm(list=ls())
graphics.off()

# Specify paths ----------------------------------------------------------------

source("code/specify_paths.R")
source("code/fn-Isq.R")
source("code/fn-uvmr.R")

# Load source data info --------------------------------------------------------

gwas <- data.table::fread("raw/gwas.csv",
                          stringsAsFactors = FALSE,
                          data.table = FALSE)

# Load instruments -------------------------------------------------------------

instruments <- data.table::fread("data/instruments_all.txt", data.table = FALSE)

# Record traits without instruments --------------------------------------------

paste0("The following traits have no instruments: ",paste(setdiff(gwas$trait, instruments$exposure), collapse = ", "))

# Record feature traits --------------------------------------------------------

features <- intersect(instruments$exposure, gwas$trait)

# Perform UVMR -----------------------------------------------------------------

results <- data.table::fread("output/results.csv", data.table = FALSE) # NULL
plei <- data.table::fread("output/plei.csv", data.table = FALSE) # NULL

for (i in features) {
  
  for (j in c("t2d","pad","cad")) {
    
    if (i!=j) {
      
      # Feature > outcome effects ----------------------------------------------
      
      tmp <- uvmr(i, j)
      results <- rbind(results, tmp[[1]])
      plei <- rbind(plei, tmp[[2]])
      
      # Outcome > feature effects ----------------------------------------------
      
      tmp <- uvmr(j, i)
      results <- rbind(results, tmp[[1]])
      plei <- rbind(plei, tmp[[2]])
      
    }
    
  }
  
}

results[,c("id.exposure","id.outcome")] <- NULL
data.table::fwrite(results,"output/results.csv")

plei[,c("id.exposure","id.outcome")] <- NULL
data.table::fwrite(plei,"output/plei.csv")