extract_snps <- function(gwas,ieugwas,snplist,label) {
  
  if (ieugwas=="") {
    
    df <- data.table::fread(paste0("data/gwas-",gwas,".txt"),
                             stringsAsFactors = FALSE,
                             data.table = FALSE)
    
    if (sum(grepl("rs",df$SNP))==0) {
      
      df$chrpos <- paste(df$chr, df$pos, sep = ":")
      
      df <- df[df$chrpos %in% snplist$chrpos,]
      
      df$SNP <- NULL
      
      df <- merge(df, snplist[,c("SNP","chrpos")], by = "chrpos")
      
    } else {
      
      df <- df[df$SNP %in% snplist$SNP,]
      
    }
    
  } else {
    
    df <- TwoSampleMR::extract_outcome_data(outcomes = ieugwas,
                                             snps = snplist$SNP,
                                             proxies = FALSE)
    
    colnames(df) <- gsub(".outcome","",colnames(df))
    
  }
  
  df <- df[,c("SNP","beta","se","pval","effect_allele","other_allele")]
  colnames(df) <- c("SNP",paste0("beta.",label),paste0("se.",label),paste0("pval.",label),paste0("effect_allele.",label),paste0("other_allele.",label))
  
  return(df)
  
}
