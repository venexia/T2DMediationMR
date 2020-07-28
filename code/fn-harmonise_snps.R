harmonise_snps <- function(snplist,t2d_instrument,t2d,feature,outcome) {
  
  ## Combine all exposures and outcome in a single dataset -------------------
  
  dat <- merge(t2d, feature, by = "SNP")
  dat <- merge(dat, outcome, by = "SNP")
  
  # Clump instrument -----------------------------------------------------------
  
  if (nrow(snplist[snplist$instrument=="feature",]) < nrow(snplist[snplist$instrument=="t2d",])) {
    
    clumped <- dat[,c("SNP","pval.feature")]

  } else {
    
    clumped <- dat[,c("SNP","pval.t2d")]

  }
  
  colnames(clumped) <- c("rsid","pval")
  
  clumped <- ieugwasr::ld_clump(clumped,
                                clump_kb = 10000,
                                clump_r2 = 0.001)
  
  dat <- dat[dat$SNP %in% clumped$rsid,]
  
  colnames(dat) <- gsub("beta.","orig.beta.",colnames(dat))
  colnames(dat) <- gsub("effect_allele.","orig.effect_allele.",colnames(dat))
  colnames(dat) <- gsub("other_allele.","orig.other_allele.",colnames(dat))
  
  ## Make all effects t2d increasing -----------------------------------------
  
  dat$beta.t2d <- ifelse(sign(dat$orig.beta.t2d)==1,dat$orig.beta.t2d,-1*dat$orig.beta.t2d)
  dat$effect_allele <- ifelse(sign(dat$orig.beta.t2d)==1,dat$orig.effect_allele.t2d,dat$orig.other_allele.t2d)
  dat$other_allele <- ifelse(sign(dat$orig.beta.t2d)==1,dat$orig.other_allele.t2d,dat$orig.effect_allele.t2d)
  
  dat$beta.feature <- ifelse(dat$effect_allele==dat$orig.effect_allele.feature & dat$other_allele==dat$orig.other_allele.feature,
                              dat$orig.beta.feature,
                              -1*dat$orig.beta.feature)
  
  dat$beta.outcome <- ifelse(dat$effect_allele==dat$orig.effect_allele.outcome & dat$other_allele==dat$orig.other_allele.outcome,
                             dat$orig.beta.outcome,
                             -1*dat$orig.beta.outcome)
  
  dat <- dat[,c("SNP","effect_allele","other_allele","beta.t2d","se.t2d","beta.feature","se.feature","beta.outcome","se.outcome")]
  
  return(dat)
  
}