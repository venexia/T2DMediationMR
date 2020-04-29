rm(list=setdiff(ls(), keep))
graphics.off()

library(magrittr)

# Load features data -----------------------------------------------------------

features <- data.table::fread("data/features_sources.csv",
                              stringsAsFactors = FALSE,
                              data.table = FALSE)

features <- features[features$filename!="",]

t2d_ins <- data.table::fread("data/instrument-t2d.txt",
                             stringsAsFactors = FALSE,
                             data.table = FALSE)

t2d_ins$neale1 <- paste(t2d_ins$chr,t2d_ins$pos,t2d_ins$effect_allele,t2d_ins$other_allele,sep = ":")
t2d_ins$neale2 <- paste(t2d_ins$chr,t2d_ins$pos,t2d_ins$other_allele,t2d_ins$effect_allele,sep = ":")
t2d_ins$ligthart1 <- paste(t2d_ins$chr,t2d_ins$pos,tolower(t2d_ins$effect_allele),tolower(t2d_ins$other_allele),sep = "_")
t2d_ins$ligthart2 <- paste(t2d_ins$chr,t2d_ins$pos,tolower(t2d_ins$other_allele),tolower(t2d_ins$effect_allele),sep = "_")

# Format each feature in turn --------------------------------------------------

# Leptin
#  Leptin adj bmi
# HOMA_IR
# HOMA_B
# albumin
# basophil

for (i in 38:nrow(features)) {
  
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
  
  # Ensure alleles are capitalized ---------------------------------------------
  
  if ("effect_allele" %in% colnames(tmp)) {
  tmp$effect_allele <- toupper(tmp$effect_allele)
  }
  
  if ("other_allele" %in% colnames(tmp)) {
  tmp$other_allele <- toupper(tmp$other_allele)
  }
  
  # Extract genome wide significant hits ---------------------------------------
  
  gws <- tmp[tmp$pval<5e-8,]
  gws$instrument <- features$trait[i]
  
  # Update chr and pos from SNP if rsIDs, chr and pos missing ------------------
  
  if (features$rsid[i]==FALSE) {
    
    if (features$consortium[i]=="Neale") {
      tmp <- tmp[tmp$SNP %in% c(t2d_ins$neale1,t2d_ins$neale2),]
    }
    if (features$consortium[i]=="Ligthart") {
      tmp <- tmp[tmp$SNP %in% c(t2d_ins$ligthart1,t2d_ins$ligthart2),]
    }
  
  } else{
    tmp <- tmp[tmp$SNP %in% t2d_ins$SNP,]
  }
  
  # Join genome-wide significant and T2D ---------------------------------------
  
  tmp$instrument <- "type 2 diabetes"
  df <- rbind(tmp,gws)
  
  # Save .txt files ------------------------------------------------------------
  
  data.table::fwrite(df,paste0(path_features_final,features$trait[i],".txt"))
  
}
