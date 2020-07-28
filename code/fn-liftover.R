liftover <- function(trait) {

  # Load input -----------------------------------------------------------------
  
  liftover_in <- data.table::fread(paste0("liftover/liftover_input-",trait,".txt"),
                                   header = FALSE,
                                   col.names = c("id36"),
                                   stringsAsFactors = FALSE,
                                   data.table = FALSE)
  # Load output ----------------------------------------------------------------
  
  liftover_out <- data.table::fread(paste0("liftover/liftover_output-",trait,".bed"),
                                    header = FALSE,
                                    col.names = c("id37"),
                                    stringsAsFactors = FALSE,
                                    data.table = FALSE)
  
  # Load error file ------------------------------------------------------------
  
  liftover_err <- data.table::fread(paste0("liftover/liftover_error-",trait,".txt"),
                                    header = FALSE,
                                    stringsAsFactors = FALSE,
                                    data.table = FALSE)
  
  # Match error message to input -----------------------------------------------
  
  liftover_err$V2 <- seq(1,2,1)
  liftover_err$V2 <- ifelse(liftover_err$V2==1,"error","id36")
  liftover_err$V3 <- rep(seq(1,nrow(liftover_err)/2,1),each = 2)
  liftover_err <- tidyr::pivot_wider(liftover_err, names_from = "V2", values_from = "V1")
  liftover_err$V3 <- NULL
  
  # Remove IDs with error from input -------------------------------------------
  
  liftover_in <- data.frame(liftover_in[!(liftover_in$id36 %in% liftover_err$id36),],
                            stringsAsFactors = FALSE)
  colnames(liftover_in) <- c("id36")
  
  # Check number of rows match -------------------------------------------------
  
  if (nrow(liftover_in)!=nrow(liftover_out)) {
    stop("Liftover input does not contain the same number of rows as liftover output!")
  }
  
  # Bind new and old IDs -------------------------------------------------------
  
  df <- cbind(liftover_in,liftover_out)
  
  # Reformat -------------------------------------------------------------------
  
  df$chr36 <- gsub(":.*","",gsub("chr","",df$id36))
  df$pos36 <- gsub(".*-","",df$id36)
  
  df$chr <- gsub(":.*","",gsub("chr","",df$id37))
  df$pos <- gsub(".*-","",df$id37)
  
  df <- df[,c("chr36","pos36","chr","pos")]
  
  # Return map -----------------------------------------------------------------
  
  return(df)
  
}