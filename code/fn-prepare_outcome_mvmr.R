prepare_outcome_mvmr <- function(outcome,snps) {
  
  if (outcome=="cad") {
    
    df <- TwoSampleMR::extract_outcome_data(snps = snps,
                                            outcomes = "ieu-a-7",
                                            proxies = FALSE)
    
    df$outcome <- j
    colnames(df) <- gsub(".outcome","",colnames(df))
    df[,c("id","originalname","outcome.deprecated","mr_keep","data_source")] <- NULL
    
  } else {
    
    df <- data.table::fread(paste0("data/gwas-",outcome,".txt"),
                            stringsAsFactors = FALSE,
                            data.table = FALSE)
    
    df <- df[df$SNP %in% snps,]
    
  }
  
  return(df)
  
}