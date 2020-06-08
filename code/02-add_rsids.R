rm(list=ls())
graphics.off()

# Specify paths ----------------------------------------------------------------

source("code/specify_paths.R")

# Load feature data ------------------------------------------------------------

gwas <- data.table::fread("raw/gwas.csv")

# Restrict to SNPs without rsIDs -----------------------------------------------

gwas <- gwas[(gwas$ieugwas=="" & 
                (gwas$consortium=="Neale") | gwas$trait %in% c("crp","pad","t2d")),]

# Extract genome wide significant hits for each GWAS ---------------------------

for (i in 1:nrow(gwas)) {
  
  df <- data.table::fread(paste0("data/gwas-",gwas$trait[i],".txt"),
                          stringsAsFactors = FALSE,
                          data.table = FALSE)
  
  df <- df[df$pval < 5e-8 & !is.na(df$pval),]
  
  if (nrow(df)>0) {
    
    # Restrict to SNPs not on X and Y chromosomes ------------------------------
    
    df <- df[df$chr %in% 1:22,]
    
    # Remove insertions --------------------------------------------------------
    
    df <- df[nchar(df$effect_allele)==1,]
    df <- df[nchar(df$other_allele)==1,]
    
    # Prepare SNPNexus input ---------------------------------------------------
    
    df$type <- "chromosome"
    df$strand <- 1
    df <- df[,c("type","chr","pos","effect_allele","other_allele","strand")]
    df <- unique(df)
    
    tmp <- df
    colnames(tmp) <- c("type","chr","pos","other_allele","effect_allele","strand")
    df <- rbind(df,tmp)
    
    if (nrow(df)<=100000) {
    
    data.table::fwrite(df,
                       paste0("data/snpnexus_input-",gwas$trait[i],".txt"),
                       sep = "\t",
                       col.names = FALSE)
      
    } else {
     
      n <- ceiling(nrow(df)/100000)
      
      for (j in 1:n) {
        
        start <- (j-1)*100000
        end <- min(j*100000,nrow(df))
        
        data.table::fwrite(df[start:end,],
                           paste0("data/snpnexus_input-",gwas$trait[i],"_",j,".txt"),
                           sep = "\t",
                           col.names = FALSE)
      }
       
    }
    
  }
  
}

# Submit files to SNPNexus to annotate with rsIDs from human assembly GRCh37
