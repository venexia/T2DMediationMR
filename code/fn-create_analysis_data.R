create_analysis_data <- function(data_t2d,data_feature,data_outcome) {
  
  ## Clump SNP list --------------------------------------------------------------
  
  snplist <- intersect(data_feature$SNP,data_t2d$SNP)
  snplist <- intersect(data_outcome$SNP,snplist)
  snplist <- unique(snplist)
  clumped <- data_feature[data_feature$SNP %in% snplist,c("SNP","pval")]
  colnames(clumped) <- c("rsid","pval")
  clumped <- ieugwasr::ld_clump(dat = clumped, clump_kb = 10000, clump_r2 = 0.001)
  
  ## Create analysis dataset -----------------------------------------------------
  
  df <- data.frame(SNP = clumped$rsid,
                   stringsAsFactors = FALSE)
  
  df <- merge(df,data_t2d[,c("SNP","effect_allele","other_allele","beta","se")],all.x = TRUE,by = "SNP")
  colnames(df) <- c("SNP","effect_allele","other_allele","beta.t2d","se.t2d")
  
  df <- merge(df,data_feature[,c("SNP","effect_allele","other_allele","beta","se")],all.x = TRUE,by = "SNP")
  colnames(df) <- c("SNP","effect_allele","other_allele","beta.t2d","se.t2d",
                    "effect_allele.feature","other_allele.feature","beta.feature.orig","se.feature")
  
  df <- merge(df,data_outcome[,c("SNP","effect_allele","other_allele","beta","se")],all.x = TRUE,by = "SNP")
  colnames(df) <- c("SNP","effect_allele","other_allele","beta.t2d","se.t2d",
                    "effect_allele.feature","other_allele.feature","beta.feature.orig","se.feature",
                    "effect_allele.outcome","other_allele.outcome","beta.outcome.orig","se.outcome")
  
  ## Harmonise direction of effect -----------------------------------------------
  
  df$beta.feature <- NA
  df$beta.feature <- ifelse(df$effect_allele==df$effect_allele.feature & df$other_allele==df$other_allele.feature,df$beta.feature.orig,df$beta.feature)
  df$beta.feature <- ifelse(df$effect_allele==df$other_allele.feature & df$other_allele==df$effect_allele.feature,-1*df$beta.feature.orig,df$beta.feature)
  
  df$beta.outcome <- NA
  df$beta.outcome <- ifelse(df$effect_allele==df$effect_allele.outcome & df$other_allele==df$other_allele.outcome,df$beta.outcome.orig,df$beta.outcome)
  df$beta.outcome <- ifelse(df$effect_allele==df$other_allele.outcome & df$other_allele==df$effect_allele.outcome,-1*df$beta.outcome.orig,df$beta.outcome)
  
  df <- df[,c("SNP","effect_allele","other_allele","beta.t2d","se.t2d","beta.feature","se.feature","beta.outcome","se.outcome")]
  
  df <- na.omit(df)
  
  return(df)
  
}