rm(list=ls())
graphics.off()

# Specify paths ----------------------------------------------------------------

source("code/specify_paths.R", echo = TRUE)

# Source functions -------------------------------------------------------------

source("code/fn-liftover.R", echo = TRUE)

# Load feature data ------------------------------------------------------------

gwas <- data.table::fread("raw/gwas.csv",
                          stringsAsFactors = FALSE,
                          data.table = FALSE)

# Restrict to GWAS using NCBI36/hg18 -------------------------------------------

gwas <- gwas[gwas$build=="NCBI36/hg18",]

# Prepare LOLIPOP consortium GWAS liftover input -------------------------------

for (i in gwas$trait) {
  
  if (gwas[gwas$trait==i,]$consortium=="LOLIPOP") {
    
    tmp <- data.table::fread(paste0("raw/gwas-",i,".txt.gz"),
                             data.table = FALSE,
                             stringsAsFactors = FALSE)
    
    tmp <- tmp[,c("Chr36","Position36")]
    
  }
  
  if (i=="isi") {
    
    tmp <- data.table::fread(paste0("raw/gwas-",i,".txt"),
                             data.table = FALSE,
                             stringsAsFactors = FALSE)
    
    tmp <- tmp[,c("chr_build36","pos_build36")]
    
  }
  
  colnames(tmp) <- c("chr","end")
  
  tmp <- tmp[tmp$end>0,]
  
  tmp$start <- tmp$end - 1
  
  tmp <- na.omit(tmp)
  
  tmp$input <- paste0("chr",tmp$chr,":",tmp$start,"-",tmp$end)
  
  tmp[,c("chr","start","end")] <- NULL
  
  data.table::fwrite(tmp,paste0("liftover/input-",i,".txt"), col.names = FALSE)
  
}