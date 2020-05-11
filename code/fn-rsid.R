rsid <- function(dat,chr,pos) {
  
  dat <- data.table::data.table(dat)
  
  # Check 'SNP' does not exist in dat ----------------------------------------
  
  if ("SNP" %in% colnames(dat)) {
    warning("The variable SNP has been overwritten.")
  }
  
  dat$SNP <- NULL
  
  # Remove X and Y chromosome SNPs -------------------------------------------
  
  dat <- dat[!(dat$chr %in% c("X","Y","x","y")),]
  
  # Ensure chromosome and position are integers in dat -----------------------
  
  dat[,c(chr)] <- as.integer(dat[,c(chr)])
  dat[,c(pos)] <- as.integer(dat[,c(pos)])
  dat$chrpos <- paste0(dat[,c(chr)],"_",dat[,c(pos)])
  
  # Loop over each chromosome --------------------------------------------------
  
  custom_link <- NULL
  
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
    
    colnames(link) <- c(chr,pos,"SNP")
    link$chrpos <- paste0(link$chr,"_",link$pos)
    
    # Merge link and dat -------------------------------------------------------
    
    link <- link[link$chrpos %in% dat$chrpos,]
    custom_link <- rbind(custom_link,link)
    
  }
  
  dat <- merge(dat,custom_link,by = c("chrpos","chr","pos"))
  dat$chrpos <- NULL
  
  return(dat)
  
}