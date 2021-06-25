rm(list=ls())
graphics.off()

# Load libraries ---------------------------------------------------------------

library(magrittr)

# Specify paths ----------------------------------------------------------------

source("code/specify_paths.R", echo = TRUE)

# Empty data frame for outcome instruments -------------------------------------

instruments <- NULL

# Format PAD GWAS --------------------------------------------------------------

df <- data.table::fread("raw/gwas-pad.csv.gz",
                        data.table = FALSE,
                        stringsAsFactors = FALSE)

df$exposure <- "pad"

df <- df[,c("ID","CHROM","POS",
            "EFFECT_ALLELE","OTHER_ALLELE","EAF",
            "BETA","SE","PVAL",
            "CASES","CONTROLS","N",
            "exposure")]

colnames(df) <- c("SNP","chr","pos",
                  "effect_allele","other_allele","eaf",
                  "beta","se","pval",
                  "ncases","ncontrols","samplesize",
                  "exposure")

data.table::fwrite(df,"data/gwas-pad.csv.gz")

instruments <- rbind(instruments,df[df$pval<5e-8,])

# Format CAD GWAS --------------------------------------------------------------

df <- data.table::fread("raw/gwas-cad.csv.gz",
                        data.table = FALSE,
                        stringsAsFactors = FALSE)

df$exposure <- "cad"
df$ncases <- 60801
df$ncontrols <- 123504
df$samplesize <- df$ncases + df$ncontrols

df <- df[,c("markername","chr","bp_hg19",
            "effect_allele","noneffect_allele","effect_allele_freq",
            "beta","se_dgc","p_dgc",
            "ncases","ncontrols","samplesize",
            "exposure")]

colnames(df) <- c("SNP","chr","pos",
                  "effect_allele","other_allele","eaf",
                  "beta","se","pval",
                  "ncases","ncontrols","samplesize",
                  "exposure")

data.table::fwrite(df,"data/gwas-cad.csv.gz")

instruments <- rbind(instruments,df[df$pval<5e-8,])


# Format linear T2D GWAS -------------------------------------------------------

df <- data.table::fread("raw/gwas-t2d_linear.csv.gz",
                        data.table = FALSE,
                        stringsAsFactors = FALSE)

df$exposure <- "t2d_linear"
df$ncase <- 24884
df$ncont <- 437996
df$samplesize <- df$ncase + df$ncont
df$mu <- df$ncase / (df$ncase + df$ncont)

df$BETA_convert <- df$BETA / (df$mu * (1 - df$mu))
df$SE_convert <- df$SE / (df$mu * (1 - df$mu))

df <- df[,c("SNP","CHR","BP",
            "ALLELE1","ALLELE0","A1FREQ",
            "BETA_convert","SE_convert","P_LINREG",
            "ncase","ncont","samplesize",
            "exposure")]

colnames(df) <- c("SNP","chr","pos",
                  "effect_allele","other_allele","eaf",
                  "beta","se","pval",
                  "ncases","ncontrols","samplesize",
                  "exposure")

data.table::fwrite(df,"data/gwas-t2d_linear.csv.gz")

instruments <- rbind(instruments,df[df$pval<5e-8,])

# Format T2D GWAS --------------------------------------------------------------

df <- data.table::fread("raw/gwas-t2d.csv.gz",
                        data.table = FALSE,
                        stringsAsFactors = FALSE)

df$exposure <- "t2d"
df$Pvalue <- as.numeric(df$Pvalue)
df$ncases <- 74124
df$ncontrols <- 824006 
df$samplesize <- df$ncases + df$ncontrols
df$SNP <- NULL

df$ref_eligible <- nchar(df$EA)==1 & nchar(df$NEA)==1 & df$EAF>0.01 & df$EAF<0.99

ref <- data.table::fread("raw/data_maf0.01_rs.bim",
                         data.table = FALSE,
                         header = FALSE, 
                         select = c(1,2,4), col.names = c("Chr","SNP","Pos"))

# Note: 
# 33% coverage in ref file of all SNPs [#SNPs = 23465132]
# 98% coverage in ref file of eligible SNPs (biallelic, MAF > 0.01) [#SNPs = 7849073]

df <- merge(df, ref, by = c("Chr","Pos"))

df <- df[,c("SNP","Chr","Pos",
            "EA","NEA","EAF",
            "Beta","SE","Pvalue",
            "ncases","ncontrols","samplesize",
            "exposure")]

colnames(df) <- c("SNP","chr","pos",
                  "effect_allele","other_allele","eaf",
                  "beta","se","pval",
                  "ncases","ncontrols","samplesize",
                  "exposure")

data.table::fwrite(df,"data/gwas-t2d.csv.gz")

instruments <- rbind(instruments,df[df$pval<5e-8,])

# Clump instruments ------------------------------------------------------------

instruments <- TwoSampleMR::format_data(instruments,
                                        type = "exposure",
                                        phenotype_col = "exposure",
                                        chr_col = "chr",
                                        pos_col = "pos",
                                        effect_allele_col = "effect_allele",
                                        other_allele_col = "other_allele",
                                        eaf_col = "eaf",
                                        beta_col = "beta",
                                        se_col = "se",
                                        pval_col = "pval",
                                        ncase_col = "ncases",
                                        ncontrol_col = "ncontrols",
                                        samplesize_col = "samplesize")

instruments <- TwoSampleMR::clump_data(instruments,
                                       clump_p1 = 1,
                                       clump_p2 = 1,
                                       clump_r2 = 0.001,
                                       clump_kb = 10000)

colnames(instruments) <- gsub(".exposure","",colnames(instruments))

instruments[,c("mr_keep","pval_origin","id")] <- NULL

# Save outcome instruments -----------------------------------------------------

data.table::fwrite(instruments,"data/instruments-outcomes.csv")
