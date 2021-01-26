rm(list=ls())

# Load libraries ---------------------------------------------------------------

library(magrittr)

# Specify paths ----------------------------------------------------------------

source("code/specify_paths.R")

# Prepare SNPNexus file --------------------------------------------------------

snpnexus <- data.table::fread("raw/snpnexus_output.txt", data.table = FALSE)
snpnexus <- snpnexus[,c("dbSNP","Chromosome","Position","REF Allele","ALT Allele (IUPAC)")]
colnames(snpnexus) <- c("SNP","chr","pos","ref_allele","alt_allele")
snpnexus <- snpnexus[snpnexus$SNP!="None",]

snpnexus$alleles <- ifelse(snpnexus$ref_allele < snpnexus$alt_allele, 
                           paste(snpnexus$ref_allele, snpnexus$alt_allele, sep = ";"),
                           paste(snpnexus$alt_allele, snpnexus$ref_allele, sep = ";"))

# Load gwas information --------------------------------------------------------

gwas <- data.table::fread("raw/gwas.csv",
                          stringsAsFactors = FALSE,
                          data.table = FALSE)

gwas <- gwas[gwas$source!="exclude_feature",]

# Extract genome wide significant hits for each GWAS ---------------------------

df <- NULL

for (i in 1:nrow(gwas)) {
  
  if (is.na(gwas$ieugwas[i])) {
    
    tmp <- data.table::fread(paste0("data/gwas-",gwas$trait[i],".txt"),
                             data.table = FALSE)
    
    tmp <- tmp[tmp$pval < 5e-8 & !is.na(tmp$pval),]
    
    tmp$SNP <- NULL
    
    tmp$effect_allele <- toupper(tmp$effect_allele)
    tmp$other_allele <- toupper(tmp$other_allele)
    
    tmp$alleles <- ifelse(tmp$effect_allele < tmp$other_allele, 
                          paste(tmp$effect_allele, tmp$other_allele, sep = ";"),
                          paste(tmp$other_allele, tmp$effect_allele, sep = ";"))
    
    tmp <- merge(tmp,snpnexus,by = c("chr","pos","alleles"))
    
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
    
    tmp$exposure <- gwas$trait[i]
    colnames(tmp) <- gsub(".exposure","",colnames(tmp))
    tmp[,c("mr_keep","pval_origin","data_source","id")] <- NULL
    
    df <- plyr::rbind.fill(df,tmp)
    
  } else {
    
    tmp <- TwoSampleMR::extract_instruments(outcomes = gwas$ieugwas[i],
                                            p1 = 5e-8,
                                            clump = FALSE)
      
    if (!is.null(tmp)) {
    
      tmp$exposure <- gwas$trait[i]
      colnames(tmp) <- gsub(".exposure","",colnames(tmp))
      tmp$ieugwas <- tmp$id
      tmp[,c("mr_keep","pval_origin","data_source","id")] <- NULL
      

      df <- plyr::rbind.fill(df,tmp)
  
    }
    
  }
  
}

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