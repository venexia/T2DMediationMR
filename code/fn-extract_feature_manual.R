extract_feature_manual <- function(trait) {
  
  df <- data.table::fread(paste0(path_features_final,trait,".txt"),
                          stringsAsFactors = FALSE,
                          data.table = FALSE)
  
  tophits <- df[df$pval < 5e-8,]
  tophits$rsid <- tophits$SNP
  tophits <- ieugwasr::ld_clump(tophits)
  
  if (nrow(tophits)>0) {
    tophits <- tophits[,intersect(colnames(tophits),c("SNP","chr","pos","effect_allele","other_allele","eaf","beta","se","pval"))]
    tophits$instrument <- trait
  }

  t2dhits <- df[df$SNP %in%  t2d_ins$SNP,]

  if (nrow(t2dhits)>0) {
    t2dhits <- t2dhits[,intersect(colnames(tophits),c("SNP","chr","pos","effect_allele","other_allele","eaf","beta","se","pval"))]
    t2dhits$instrument <- "type 2 diabetes"
  }
  
  df <- rbind(tophits,t2dhits)
  
  df$trait <- trait
  
  return(df)
  
}