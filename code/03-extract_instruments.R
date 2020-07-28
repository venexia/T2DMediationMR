rm(list=ls())
graphics.off()

# Load libraries ---------------------------------------------------------------

library(magrittr)

# Specify paths ----------------------------------------------------------------

source("code/specify_paths.R")

# Load feature data ------------------------------------------------------------

gwas <- data.table::fread("raw/gwas.csv",
                          stringsAsFactors = FALSE,
                          data.table = FALSE)

snpnexus <- gwas[gwas$rsid==FALSE,]$trait

# Extract genome wide significant hits for each GWAS ---------------------------

df <- NULL

for (i in 1:nrow(gwas)) {
  
  if (gwas$ieugwas[i]=="") {
    
    tmp <- data.table::fread(paste0("data/gwas-",gwas$trait[i],".txt"),
                             stringsAsFactors = FALSE,
                             data.table = FALSE)
  
    tmp <- tmp[tmp$pval < 5e-8 & !is.na(tmp$pval),]
    
    tmp$effect_allele <- toupper(tmp$effect_allele)
    tmp$other_allele <- toupper(tmp$other_allele)
    
    if (gwas$trait[i] %in% snpnexus) {
      
      snpnexus_files <- list.files(path = "snpnexus/", pattern = paste0("snpnexus_output-",gwas$trait[i]))
      
      snpnexus_output <- NULL
      
      for (j in snpnexus_files) {
        
        tmp_snpnexus <- data.table::fread(paste0("snpnexus/",j),
                                          stringsAsFactors = FALSE,
                                          data.table = FALSE)
        
        snpnexus_output <- rbind(snpnexus_output,tmp_snpnexus)
        
      }
      
      snpnexus_output <-  snpnexus_output[,c("dbSNP","Chromosome","Position","REF Allele","ALT Allele (IUPAC)")]
      colnames(snpnexus_output) <- c("SNP","chr","pos","ref_allele","alt_allele")
      snpnexus_output <- snpnexus_output[snpnexus_output$SNP!="None",]
      
      snpnexus_output$alleles <- ifelse(snpnexus_output$ref_allele < snpnexus_output$alt_allele, 
                                        paste(snpnexus_output$ref_allele, snpnexus_output$alt_allele, sep = ";"),
                                        paste(snpnexus_output$alt_allele, snpnexus_output$ref_allele, sep = ";"))
      
      tmp$SNP <- NULL
      
      tmp$alleles <- ifelse(tmp$effect_allele < tmp$other_allele, 
                            paste(tmp$effect_allele, tmp$other_allele, sep = ";"),
                            paste(tmp$other_allele, tmp$effect_allele, sep = ";"))
      
      tmp <- merge(tmp,snpnexus_output,by = c("chr","pos","alleles"))
      
    }
    
    if (nrow(tmp)>0) {
      
      # When rsID is in data twice, keep lowest p-value ------------------------
      
      tmp <- tmp %>%
        dplyr::group_by(SNP) %>%
        dplyr::slice(which.min(pval))
      
      tmp <- suppressWarnings(TwoSampleMR::format_data(dat = tmp,
                                                       type = "exposure",
                                                       snps = NULL,
                                                       header = TRUE,
                                                       phenotype_col = "exposure",
                                                       snp_col = "SNP",
                                                       beta_col = "beta",
                                                       se_col = "se",
                                                       eaf_col = "eaf",
                                                       effect_allele_col = "effect_allele",
                                                       other_allele_col = "other_allele",
                                                       pval_col = "pval",
                                                       samplesize_col = "samplesize",
                                                       chr_col = "chr",
                                                       pos_col = "pos",
                                                       log_pval = FALSE))
      
      # tmp <- TwoSampleMR::clump_data(tmp,
      #                                clump_kb = 10000,
      #                                clump_r2 = 0.001,
      #                                clump_p1 = 5e-8,
      #                                clump_p2 = 5e-8,
      #                                pop = "EUR")
      
      tmp$exposure <- gwas$trait[i]
      colnames(tmp) <- gsub(".exposure","",colnames(tmp))
      tmp[,c("mr_keep","pval_origin","data_source","id")] <- NULL
      
      # data.table::fwrite(tmp,paste0("data/instrument-",gwas$trait[i],".txt"))
      
      df <- plyr::rbind.fill(df,tmp)
      
    }
    
  } else {
    
    tmp <- suppressWarnings(TwoSampleMR::extract_instruments(outcomes = gwas$ieugwas[i], 
                                                             p1 = 5e-8,
                                                             clump = FALSE))
    
    # tmp <- suppressWarnings(TwoSampleMR::extract_instruments(outcomes = gwas$ieugwas[i], 
    #                                                          p1 = 5e-8,
    #                                                          clump = TRUE,
    #                                                          p2 = 5e-8,
    #                                                          r2 = 0.001,
    #                                                          kb = 10000))
    
    if (!is.null(tmp)) {
      
      tmp$exposure <- gwas$trait[i]
      colnames(tmp) <- gsub(".exposure","",colnames(tmp))
      tmp$ieugwas <- tmp$id
      tmp[,c("mr_keep","pval_origin","data_source","id")] <- NULL
      
      # data.table::fwrite(tmp,paste0("data/instrument-",gwas$trait[i],".txt"))
      
      df <- plyr::rbind.fill(df,tmp)
      
    }
  }
  
}

# Add Udler type 2 diabetes instrument -----------------------------------------

udler <- data.table::fread("raw/T2D_259snps_RISK_ALLELE.txt",
                           stringsAsFactors = FALSE,
                           data.table = FALSE)

tmp <- data.table::fread("data/gwas-t2d.txt",
                         stringsAsFactors = FALSE,
                         data.table = FALSE)

udler <- tidyr::separate(data = udler,
                         col = VAR_ID_hg19,
                         into = c("chr","pos","ref_allele","alt_allele"),
                         sep = "_",
                         remove = FALSE)

udler$chrpos <- paste(udler$chr,udler$pos,sep = ":")

tmp$chrpos <- tmp$SNP

tmp$SNP <- NULL

tmp <- tmp[tmp$chrpos %in% udler$chrpos,]

tmp <- merge(tmp,udler[,c("chrpos","SNP")])

tmp$chrpos <- NULL

tmp$exposure <- "t2d_udler"

df <- plyr::rbind.fill(df,tmp)

# Generate chromosome position annotation for instruments ----------------------

df$chrpos <- paste(df$chr, df$pos, sep = ":")

# Make all effects positive ---------------------------------------------------

df$beta.original <- df$beta
df$eaf.original <- df$eaf
df$effect_allele.original <- df$effect_allele
df$other_allele.original <- df$other_allele
df[,c("beta","eaf","effect_allele","other_allele")] <- NULL

df <- df %>%
  dplyr::mutate(beta = ifelse(sign(beta.original)==-1, -1*beta.original, beta.original)) %>%
  dplyr::mutate(effect_allele = ifelse(sign(beta.original)==-1, other_allele.original, effect_allele.original)) %>%
  dplyr::mutate(other_allele = ifelse(sign(beta.original)==-1, effect_allele.original, other_allele.original)) %>%
  dplyr::mutate(eaf = ifelse(beta.original <= 0, 1-eaf.original, eaf.original))

df[,c("beta.original","eaf.original","effect_allele.original","other_allele.original")] <- NULL

# Save all instruments ---------------------------------------------------------

data.table::fwrite(df,"data/instruments_all.txt")