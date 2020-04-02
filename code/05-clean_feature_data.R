rm(list=setdiff(ls(), keep))
graphics.off()

# Load features data -----------------------------------------------------------

features <- data.table::fread("data/features_sources.csv",
                              stringsAsFactors = FALSE,
                              data.table = FALSE)

features <- features[features$filename!="",]

# Format each feature in turn --------------------------------------------------

for (i in 9:nrow(features)) {

  # State feature under consideration ------------------------------------------
  
  print(paste0("Feature: ",features[i,]$trait_long))
    
  # Load data ------------------------------------------------------------------
  
  tmp <- data.table::fread(paste0(path_features_processed,features$trait[i],features$ext[i]),
                           skip = features$skip[i],
                           stringsAsFactors = FALSE,
                           data.table = FALSE)
  
  # Name columns if header==FALSE and consortium=="GIANT" ----------------------
  
  if (features$header[i]==FALSE & features$consortium[i]=="GIANT") {
    colnames(tmp) <- c("MarkerName","Allele1","Allele2","FreqAllele1HapMapCEU","b","se","p","n")
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
                                                   ncase_col = ifelse(features$n_case[i]=="","ncase",features$n_case[i]),
                                                   ncontrol_col = ifelse(features$n_control[i]=="","ncontrol",features$n_control[i]),
                                                   samplesize_col = ifelse(features$n[i]=="","samplesize",features$n[i]),
                                                   z_col = ifelse(features$zscore[i]=="","z",features$zscore[i]),
                                                   chr_col = ifelse(features$chr_id[i]=="","chr",features$chr_id[i]),
                                                   pos_col = ifelse(features$chr_pos[i]=="","pos",features$chr_pos[i])))
  
  # Rename columns -------------------------------------------------------------
  
  colnames(tmp) <- gsub(pattern = ".exposure","",colnames(tmp))
  names(tmp)[names(tmp)=="exposure"] <- "trait"
  
  # Remove unnecessary columns -------------------------------------------------
  
  tmp[,c("mr_keep","pval_origin","id")] <- NULL
  
  # Save .txt files ------------------------------------------------------------------
  
  data.table::fwrite(tmp,paste0(path_features_final,features$trait[i],".txt"))
  
}
