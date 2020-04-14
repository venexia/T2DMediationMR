rm(list=setdiff(ls(), keep))
graphics.off()

# Load GWAS --------------------------------------------------------------------

df <- data.table::fread("raw/Mahajan.NatGenet2018b.t2d.European.txt",
                        stringsAsFactors = FALSE,
                        data.table = TRUE)

# Format GWAS ------------------------------------------------------------------

df <- suppressWarnings(TwoSampleMR::format_data(df,
                                                type = "exposure",
                                                snps = NULL,
                                                header = TRUE,
                                                snp_col = "SNP",
                                                beta_col = "Beta",
                                                se_col = "SE",
                                                eaf_col = "EAF",
                                                effect_allele_col = "EA",
                                                other_allele_col = "NEA",
                                                pval_col = "Pvalue",
                                                samplesize_col = "Neff",
                                                chr_col = "Chr",
                                                pos_col = "Pos"))

# Tidy data --------------------------------------------------------------------

colnames(df) <- gsub(pattern = ".exposure","",colnames(df))
df$trait <- "type 2 diabetes"
df[,c("exposure","id","mr_keep","pval_origin")] <- NULL

# Mark Udler instruments -------------------------------------------------------

ins_udler <- data.table::fread("raw/T2D_259snps_RISK_ALLELE.txt",
                               select = c("SNP","VAR_ID_hg19"))

ins_udler <- cbind(ins_udler,stringr::str_split_fixed(ins_udler$VAR_ID_hg19,"_",3))
ins_udler[,c("VAR_ID_hg19","V3")] <- NULL
colnames(ins_udler) <- c("SNP","chr","pos")

ins_udler$chr <- as.numeric(ins_udler$chr)
ins_udler$pos <- as.numeric(ins_udler$pos)

ins_udler$rsid <- ins_udler$SNP
ins_udler$SNP <- paste0(ins_udler$chr,":",ins_udler$pos)

df$instrument_udler <- df$SNP %in% ins_udler$SNP

# Save data --------------------------------------------------------------------

data.table::fwrite(df,"data/t2d.txt", row.names = FALSE)

# Add rsIDs to GWAS ------------------------------------------------------------

df <- rsid(df,"chr","pos")

# Save data --------------------------------------------------------------------

data.table::fwrite(df,"data/t2d.txt", row.names = FALSE)

# Mark instruments -------------------------------------------------------------

ins <- df[df$pval < 5e-8 & !is.na(df$pval),]
ins$rsid <- ins$SNP
ins <- ieugwasr::ld_clump(dat = ins, clump_kb = 10000, clump_r2 = 0.001)

df$instrument <- df$SNP %in% ins$rsid

# Save data --------------------------------------------------------------------

data.table::fwrite(df,"data/t2d.txt", row.names = FALSE)
