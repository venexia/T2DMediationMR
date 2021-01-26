rm(list=ls())

annotate <- c("t2d","alp","alt","ast","ggt")

df <- NULL

for (i in annotate) {
  
  tmp <- data.table::fread(paste0("data/gwas-",i,".txt"),
                           data.table = FALSE,
                           select = c("chr","pos","effect_allele","other_allele", "pval"),
                           stringsAsFactors = FALSE)
  
  tmp <- tmp[tmp$pval<5e-8,]
  
  tmp$pval <- NULL
  
  df <- rbind(df, tmp)
  
}

# Restrict to SNPs not on X and Y chromosomes ----------------------------------

df <- df[df$chr %in% 1:22,]

# Remove indels ----------------------------------------------------------------

df <- df[nchar(df$effect_allele)==1,]
df <- df[nchar(df$other_allele)==1,]

# Prepare SNPNexus input -------------------------------------------------------

df$type <- "chromosome"
df$strand <- 1
df <- df[,c("type","chr","pos","effect_allele","other_allele","strand")]
df <- unique(df)

tmp <- df
colnames(tmp) <- c("type","chr","pos","other_allele","effect_allele","strand")
df <- rbind(df,tmp)

# Save SNPNexus input ----------------------------------------------------------

data.table::fwrite(df,"data/snpnexus_input.txt",sep = "\t",col.names = FALSE)

