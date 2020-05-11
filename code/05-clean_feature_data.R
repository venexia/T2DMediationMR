rm(list=setdiff(ls(), keep))
graphics.off()

# Load features data -----------------------------------------------------------

features <- data.table::fread("data/features_sources.csv",
                              stringsAsFactors = FALSE,
                              data.table = FALSE)

features <- features[features$ieugwas=="",]

# Format each feature in turn --------------------------------------------------

for (i in 15:nrow(features)) {
  
  # State feature under consideration ------------------------------------------
  
  print(paste0("Feature: ",features$trait_long[i]))
  
  # Load data ------------------------------------------------------------------
  
  tmp <- data.table::fread(paste0(path_features_processed,features$trait[i],features$ext[i]),
                           skip = features$skip[i],
                           stringsAsFactors = FALSE,
                           data.table = FALSE)
  
  # Name columns if header missing ---------------------------------------------
  
  if (features$header[i]==FALSE) {
    f <- features[i,c(10:19)]
    f <- f[!is.na(f) & f!=""]
    colnames(tmp) <- f
  }
  
  # Name trait -----------------------------------------------------------------
  
  tmp$Phenotype <- features$trait[i]
  
  # Format data ----------------------------------------------------------------
  
  tmp <- suppressWarnings(TwoSampleMR::format_data(tmp,
                                                   type = "exposure",
                                                   snps = NULL,
                                                   header = TRUE,
                                                   phenotype_col = "Phenotype",
                                                   snp_col = ifelse(features$snp[i]=="","SNP",features$snp[i]),
                                                   beta_col = ifelse(features$beta[i]=="","beta",features$beta[i]),
                                                   se_col = ifelse(features$se[i]=="","se",features$se[i]),
                                                   eaf_col = ifelse(features$eaf[i]=="","eaf",features$eaf[i]),
                                                   effect_allele_col = ifelse(features$effect_allele[i]=="","effect_allele",features$effect_allele[i]),
                                                   other_allele_col = ifelse(features$other_allele[i]=="","other_allele",features$other_allele[i]),
                                                   pval_col = ifelse(features$pvalue[i]=="","pval",features$pvalue[i]),
                                                   samplesize_col = ifelse(features$n[i]=="","samplesize",features$n[i]),
                                                   chr_col = ifelse(features$chr_id[i]=="","chr",features$chr_id[i]),
                                                   pos_col = ifelse(features$chr_pos[i]=="","pos",features$chr_pos[i])))
  
  # Rename columns -------------------------------------------------------------
  
  colnames(tmp) <- gsub(pattern = ".exposure","",colnames(tmp))
  names(tmp)[names(tmp)=="exposure"] <- "trait"
  
  # Remove unnecessary columns -------------------------------------------------
  
  tmp[,c("mr_keep","pval_origin","id")] <- NULL
  
  # Remove empty columns -------------------------------------------------------
  
  for (j in colnames(tmp)) {
    if (class(tmp[,c(j)])=="character") {
      if (any(tmp[,c(j)]!="")==FALSE) {
        tmp[,c(j)] <- NULL
      } 
    } else {
      if (any(!is.na(tmp[,c(j)]))==FALSE) {
        tmp[,c(j)] <- NULL
      }
    }      
  }
  
  # Update chr and pos from SNP if rsIDs, chr and pos missing ------------------
  
  if (features$rsid[i]==FALSE) {
    
    if (features$consortium[i]=="Neale") {
      tmp <- tmp %>%
        tidyr::separate(SNP, c("chr", "pos","allele1","allele2"), ":", remove = FALSE)
      tmp$variant <- tmp$SNP
      tmp$SNP <- NULL
    }
    if (features$consortium[i]=="Ligthart") {
      tmp <- tmp %>%
        tidyr::separate(SNP, c("chr", "pos","allele1","allele2"), "_", remove = FALSE)
      tmp$variant <- tmp$SNP
      tmp$SNP <- NULL
    }
    
    tmp$allele1 <- toupper(tmp$allele1)
    tmp$allele2 <- toupper(tmp$allele2)
    
    tmp$other_allele <- NA
    tmp$other_allele <- ifelse(tmp$allele1==tmp$effect_allele,tmp$allele2,tmp$other_allele)
    tmp$other_allele <- ifelse(tmp$allele2==tmp$effect_allele,tmp$allele1,tmp$other_allele)
    
    tmp[,c("allele1","allele2")] <- NULL
    
  }
  
  # Add rsIDs if missing -------------------------------------------------------
  
  if (features$rsid[i]==FALSE) {
    tmp <- rsid(tmp,"chr","pos")
  }
  
  # Ensure alleles are capitalized ---------------------------------------------
  
  tmp$effect_allele <- toupper(tmp$effect_allele)
  tmp$other_allele <- toupper(tmp$other_allele)
  tmp <- tmp[!is.na(tmp$SNP),]
  
  # Save .txt files ------------------------------------------------------------
  
  data.table::fwrite(tmp,paste0(path_features_final,features$trait[i],".txt"))
  
}