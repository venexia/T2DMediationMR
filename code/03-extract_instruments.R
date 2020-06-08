rm(list=ls())
graphics.off()

# Specify paths ----------------------------------------------------------------

source("code/specify_paths.R")

# Load feature data ------------------------------------------------------------

gwas <- data.table::fread("raw/gwas.csv",
                          stringsAsFactors = FALSE,
                          data.table = FALSE)

snpnexus <-  gwas[(gwas$ieugwas=="" & 
                     (gwas$consortium=="Neale") | gwas$trait %in% c("crp","pad","t2d")),]$trait

# Extract genome wide significant hits for each GWAS ---------------------------

df <- NULL

for (i in 1:nrow(gwas)) {
  
  if (gwas$ieugwas[i]=="") {
    
    tmp <- data.table::fread(paste0("data/gwas-",gwas$trait[i],".txt"),
                             stringsAsFactors = FALSE,
                             data.table = FALSE)
    
    
    tmp <- tmp[tmp$pval < 5e-8 & !is.na(tmp$pval),]
    
    if (gwas$trait[i] %in% snpnexus) {
      
      snpnexus_files <- list.files(path = "raw/", pattern = paste0("snpnexus_output-",gwas$trait[i]))
      
      snpnexus_output <- NULL
      
      for (j in snpnexus_files) {
        
        tmp_snpnexus <- data.table::fread(paste0("raw/",j),
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
      
      tmp1 <- tmp[,c("SNP","pval")]
      colnames(tmp1) <- c("SNP","P")
      data.table::fwrite(tmp1,paste0("data/plink_",gwas$trait[i],".assoc"),sep="\t",na = "")
      system(paste0("./plink_mac_20200428/plink --bfile raw/data_maf0.01_rs --clump data/plink_",gwas$trait[i],".assoc --clump-kb 10000 --clump-r2 0.001 --clump-p1 1 --clump-p2 1 --out data/plink_",gwas$trait[i])) 
      
      tmp1 <- data.table::fread(paste0("data/plink_",gwas$trait[i],".clumped"),
                                stringsAsFactors = FALSE,
                                data.table = FALSE)
      
      tmp <- tmp[tmp$SNP %in% tmp1$SNP,]
      
      # data.table::fwrite(tmp,paste0("data/instrument-",gwas$trait[i],".txt"))
      
      df <- plyr::rbind.fill(df,tmp)
      
    }
    
  } else {
    
    tmp <- TwoSampleMR::extract_instruments(gwas$ieugwas[i])
    
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

# Save all instruments ---------------------------------------------------------

data.table::fwrite(df,"data/instruments_all.txt")