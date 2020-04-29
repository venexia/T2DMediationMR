extract_feature <- function(trait,id) {
  
  tophits <- as.data.frame(ieugwasr::tophits(id = id))
  
  if (nrow(tophits)>0) {
    tophits <- tophits[,c("rsid","chr","position","ea","nea","eaf","beta","se","p")]
    tophits$instrument <- trait
  }

  t2dhits <- as.data.frame(ieugwasr::associations(variants = t2d_ins$SNP, id = id, proxies = FALSE))

  if (nrow(t2dhits)>0) {
    t2dhits <- t2dhits[,c("rsid","chr","position","ea","nea","eaf","beta","se","p")]
    t2dhits$instrument <- "type 2 diabetes"
  }
  
  df <- rbind(tophits,t2dhits)
  
  df$trait <- trait
  
  colnames(df) <- c("SNP","chr","pos","effect_allele","other_allele","eaf","beta","se","pval","instrument","trait")
 
  return(df)
  
}