rm(list=ls())
graphics.off()

# Specify paths ----------------------------------------------------------------

source("code/specify_paths.R", echo = TRUE)

# List GWAS with NCBI36/hg18 genome annotation ---------------------------------

update <- c("alp","alt","ast","ggt")

# Create empty liftover input file ---------------------------------------------

liftover_in <- NULL

# Prepare liftover input file by looping over eligible GWAS --------------------

for (i in update) {
  
  ## Load GWAS with NCBI36/hg18 genome annotation ------------------------------
  
  tmp <- data.table::fread(paste0("raw/gwas-",i,".txt.gz"),
                           data.table = FALSE,
                           stringsAsFactors = FALSE)
  
  ## Keep copy of GWAS to save reloading later ---------------------------------
  
  assign(i,tmp)
  
  ## Restrict columns to genome annotation -------------------------------------
  
  tmp <- tmp[,c("Chr36","Position36")]
  colnames(tmp) <- c("chr","end")
  tmp <- tmp[tmp$end>0,]
  
  ## Derive start positions ----------------------------------------------------
  
  tmp$start <- tmp$end - 1
  
  ## Format for liftover software -----------------------------------------------
  
  tmp <- na.omit(tmp)
  tmp$id36 <- paste0("chr",tmp$chr,":",tmp$start,"-",tmp$end)
  tmp[,c("chr","start","end")] <- NULL

  ## Apend to liftover input file ----------------------------------------------
  
  liftover_in <- rbind(liftover_in,tmp)
  
}

# Remove duplicates from liftover input file and save --------------------------

liftover_in <- unique(liftover_in)
data.table::fwrite(liftover_in, "data/liftover_input.txt", col.names = FALSE)

# ---------------------------------------------------------------------------------------- #
# Lift genome annotations to GRCh27/hg19 using https://genome.ucsc.edu/cgi-bin/hgLiftOver  #
# ---------------------------------------------------------------------------------------- #

# Load liftover output ---------------------------------------------------------

liftover_out <- data.table::fread("raw/liftover_output.bed",
                                  header = FALSE,
                                  col.names = c("id37"),
                                  stringsAsFactors = FALSE,
                                  data.table = FALSE)

# Load liftover errors ---------------------------------------------------------

liftover_err <- data.table::fread("raw/liftover_errors.txt",
                                  header = FALSE,
                                  stringsAsFactors = FALSE,
                                  data.table = FALSE)

# Match liftover error to liftover input ---------------------------------------

liftover_err$V2 <- seq(1,2,1)
liftover_err$V2 <- ifelse(liftover_err$V2==1,"error","id36")
liftover_err$V3 <- rep(seq(1,nrow(liftover_err)/2,1),each = 2)
liftover_err <- tidyr::pivot_wider(liftover_err, names_from = "V2", values_from = "V1")
liftover_err$V3 <- NULL

# Remove IDs with error from liftover input ------------------------------------

liftover_in <- data.frame(liftover_in[!(liftover_in$id36 %in% liftover_err$id36),],
                          stringsAsFactors = FALSE)

colnames(liftover_in) <- c("id36")

# Create liftover mapping file -------------------------------------------------

map <- cbind(liftover_in,liftover_out)

# Reformat mapping file --------------------------------------------------------

map$chr36 <- gsub(":.*","",gsub("chr","",map$id36))
map$pos36 <- gsub(".*-","",map$id36)

map$chr <- gsub(":.*","",gsub("chr","",map$id37))
map$pos <- gsub(".*-","",map$id37)

map <- map[,c("chr36","pos36","chr","pos")]

data.table::fwrite(map, "data/liver_enzyme_snps.csv")

# Apply map to each GWAS -------------------------------------------------------

for (i in update) {
  
  ## Load GWAS -----------------------------------------------------------------
  
  tmp <- get(i)
  
  ## Reformat GWAS -------------------------------------------------------------
  
  tmp$exposure <- i
  
  tmp <- tmp[,c("Chr36","Position36",
                "Effect_allele","Other_allele",
                "Effect","StdErr","P_value",
                "N",
                "exposure")]
  
  colnames(tmp) <- c("chr36","pos36",
                     "effect_allele","other_allele",
                     "beta","se","pval",
                     "samplesize",
                     "exposure")
  
  # Merge GWAS with mapping file -----------------------------------------------
  
  tmp <- merge(tmp, map, by = c("chr36","pos36"))
  tmp[,c("chr36","pos36")] <- NULL
  
  # Save GWAS with GRCh27/hg19 genome annotation -------------------------------
  
  data.table::fwrite(tmp,paste0("data/gwas-",i,".txt"))
  
}

# Format t2d GWAS --------------------------------------------------------------

tmp <- data.table::fread("raw/gwas-t2d.txt",
                         data.table = FALSE,
                         stringsAsFactors = FALSE)

tmp$exposure <- "t2d"
tmp$Pvalue <- as.numeric(tmp$Pvalue)
tmp$ncases <- 74124
tmp$ncontrols <- 824006 
tmp$samplesize <- tmp$ncases + tmp$ncontrols

tmp <- tmp[,c("SNP","Chr","Pos",
              "EA","NEA","EAF",
              "Beta","SE","Pvalue",
              "ncases","ncontrols","samplesize",
              "exposure")]

colnames(tmp) <- c("SNP","chr","pos",
                   "effect_allele","other_allele","eaf",
                   "beta","se","pval",
                   "ncases","ncontrols","samplesize",
                   "exposure")

data.table::fwrite(tmp,"data/gwas-t2d.txt")

# Format pad GWAS --------------------------------------------------------------

tmp <- data.table::fread("raw/gwas-pad.csv",
                         data.table = FALSE,
                         stringsAsFactors = FALSE)

tmp$exposure <- "pad"

tmp <- tmp[,c("ID","CHROM","POS",
              "EFFECT_ALLELE","OTHER_ALLELE","EAF",
              "BETA","SE","PVAL",
              "CASES","CONTROLS","N",
              "exposure")]

colnames(tmp) <- c("SNP","chr","pos",
                   "effect_allele","other_allele","eaf",
                   "beta","se","pval",
                   "ncases","ncontrols","samplesize",
                   "exposure")

data.table::fwrite(tmp,"data/gwas-pad.txt")