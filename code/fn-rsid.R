rsid <- function(dat,chr,pos) {
  
  dat <- data.table::data.table(dat)
  
  # Check 'SNP' does not exist in dat ----------------------------------------
  
  if ("SNP" %in% colnames(dat)) {
    warning("The variable SNP has been overwritten.")
  }
  
  dat$SNP <- NA
  
  # Ensure chromosome and position are integers in dat -----------------------
  
  dat[,c(chr)] <- as.integer(dat[,c(chr)])
  dat[,c(pos)] <- as.integer(dat[,c(pos)])
  
  # Loop over each chromosome --------------------------------------------------
  
  for (i in 1:22) {
    
    # Load link for chromosome i -----------------------------------------------
    
    link <- data.table::fread(paste0(path_project_box,"snps/snps_chr",i),
                              stringsAsFactors = FALSE,
                              data.table = TRUE)
    
    # Ensure chromosome and position are integers in link ----------------------
    
    link$`#chrom` <- gsub("chr","",link$`#chrom`)
    link$`#chrom` <- as.integer(link$`#chrom`)
    link$chromEnd <- as.integer(link$chromEnd)
    
    # Rename columns in link to match dat --------------------------------------
    
    colnames(link) <- c(chr,pos,"tmp_SNP")
    
    # Merge link and dat -------------------------------------------------------
    
    dat <- data.table::merge.data.table(dat,link,by=c(chr,pos),all.x = TRUE)
    dat$SNP <- ifelse(dat$chr==i,dat$tmp_SNP,dat$SNP)
    dat$tmp_SNP <- NULL
    
  }
  
  return(dat)
  
}