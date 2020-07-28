rm(list=ls())
graphics.off()

# Specify paths ----------------------------------------------------------------

source("code/specify_paths.R", echo = TRUE)

# Source functions -------------------------------------------------------------

source("code/fn-liftover.R", echo = TRUE)

# Load feature data ------------------------------------------------------------

gwas <- data.table::fread("raw/gwas.csv")

# Restrict to GWAS in IEU GWAS database ----------------------------------------

gwas <- gwas[gwas$ieugwas=="",]

# Format LOLIPOP consortium GWAS -----------------------------------------------

for (i in gwas[gwas$consortium=="LOLIPOP",]$trait) {
  
  tmp <- data.table::fread(paste0("raw/gwas-",i,".txt.gz"),
                           data.table = FALSE,
                           stringsAsFactors = FALSE)
  
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
  
  map <- liftover(trait = i)
  tmp <- merge(tmp, map, by = c("chr36","pos36"))
  
  tmp[,c("chr36","pos36")] <- NULL
  
  data.table::fwrite(tmp,paste0("data/gwas-",i,".txt"))
  
}

# Format Neale lab GWAS --------------------------------------------------------

for (i in gwas[gwas$consortium=="Neale",]$trait) {
  
  system(paste0("gunzip -c raw/gwas-",i,".tsv.bgz > data/gwas-",i,".txt"))
  
  tmp <- data.table::fread(paste0("data/gwas-",i,".txt"),
                           data.table = FALSE,
                           stringsAsFactors = FALSE)
  
  tmp$exposure <- i
  
  tmp <- tidyr::separate(data = tmp,
                         col = variant,
                         into = c("chr","pos","effect_allele","other_allele"),
                         sep = ":",
                         remove = FALSE)
  
  tmp <- tmp[,c("variant","chr","pos",
                "effect_allele","other_allele",
                "beta","se","pval",
                "exposure")]
  
  colnames(tmp) <- c("SNP","chr","pos",
                     "effect_allele","other_allele",
                     "beta","se","pval",
                     "exposure")
  
  data.table::fwrite(tmp,paste0("data/gwas-",i,".txt"))
  
}

# Format GIANT consortium GWAS -------------------------------------------------

for (i in gwas[gwas$consortium=="GIANT",]$trait) {
  
  tmp <- data.table::fread(paste0("raw/gwas-",i,".txt.gz"),
                           data.table = FALSE,
                           stringsAsFactors = FALSE)
  
  tmp$exposure <- i
  
  tmp <- tmp[,c("SNP","CHR","POS",
                "Tested_Allele","Other_Allele","Freq_Tested_Allele_in_HRS",
                "BETA","SE","P",
                "N",
                "exposure")]
  
  colnames(tmp) <- c("SNP","chr","pos",
                     "effect_allele","other_allele","eaf",
                     "beta","se","pval",
                     "samplesize",
                     "exposure")
  
  data.table::fwrite(tmp,paste0("data/gwas-",i,".txt"))
  
}


# Format crp GWAS --------------------------------------------------------------

tmp <- data.table::fread("raw/gwas-crp.txt",
                         data.table = FALSE,
                         stringsAsFactors = FALSE)

tmp$exposure <- "crp"

tmp <- tidyr::separate(data = tmp,
                       col = VAR_ID,
                       into = c("chr","pos","ref_allele","alt_allele"),
                       sep = "_",
                       remove = FALSE)

tmp$other_allele <- NA
tmp$other_allele <- ifelse(tmp$Effect_Allele_PH==tmp$ref_allele,tmp$alt_allele,tmp$other_allele)
tmp$other_allele <- ifelse(tmp$Effect_Allele_PH==tmp$alt_allele,tmp$ref_allele,tmp$other_allele)

tmp <- tmp[,c("VAR_ID","chr","pos",
              "Effect_Allele_PH","other_allele",
              "BETA","SE","P_VALUE",
              "exposure")]

colnames(tmp) <- c("SNP","chr","pos",
                   "effect_allele","other_allele",
                   "beta","se","pval",
                   "exposure")

data.table::fwrite(tmp,"data/gwas-crp.txt")

# Format isi GWAS --------------------------------------------------------------

tmp <- data.table::fread("raw/gwas-isi_adjbmi.txt",
                         data.table = FALSE,
                         stringsAsFactors = FALSE)

tmp$exposure <- "isi_adjbmi"

tmp <- tmp[,c("chr_build36","pos_build36",
              "effect_allele","other_allele",
              "effect","stderr","pvalue",
              "exposure")]

colnames(tmp) <- c("chr36","pos36",
                   "effect_allele","other_allele",
                   "beta","se","pval",
                   "exposure")

map <- liftover(trait = "isi_adjbmi")
tmp <- merge(tmp, map, by = c("chr36","pos36"))

tmp[,c("chr36","pos36")] <- NULL

tmp$effect_allele <- toupper(tmp$effect_allele)
tmp$other_allele <- toupper(tmp$other_allele)

data.table::fwrite(tmp,"data/gwas-isi_adjbmi.txt")

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