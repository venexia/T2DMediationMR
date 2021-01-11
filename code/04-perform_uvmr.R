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

results <- data.table::fread("output/results_201126.csv", data.table = FALSE)

plei <- data.table::fread("output/plei_201126.csv", data.table = FALSE)

for (i in features[16:71]) {
  
  for (j in c("t2d","pad","cad")) {
    
    if (i!=j) {
      
      # Feature > outcome effects ----------------------------------------------
      
      tmp <- uvmr(i, j)
      
      tmp[[1]]$id.exposure <- NULL
      tmp[[1]]$id.outcome <- NULL
      results <- rbind(results, tmp[[1]])
      
      tmp[[2]]$id.exposure <- NULL
      tmp[[2]]$id.outcome <- NULL
      plei <- rbind(plei, tmp[[2]])
      
      # Outcome > feature effects ----------------------------------------------
      
      tmp <- uvmr(j, i)
      
      tmp[[1]]$id.exposure <- NULL
      tmp[[1]]$id.outcome <- NULL
      results <- rbind(results, tmp[[1]])
      
      tmp[[2]]$id.exposure <- NULL
      tmp[[2]]$id.outcome <- NULL
      plei <- rbind(plei, tmp[[2]])
      
    }
    
  }
  
}


data.table::fwrite(results,"output/results.csv")
data.table::fwrite(plei,"output/plei.csv")