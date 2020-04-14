rm(list=setdiff(ls(), keep))
graphics.off()

# Clean CAD GWAS ---------------------------------------------------------------

cad <- data.table::fread("raw/cad.add.160614.website.txt",
                         stringsAsFactors = FALSE,
                         data.table = FALSE)

cad <- suppressWarnings(TwoSampleMR::format_data(cad,
                                                 type = "exposure",
                                                 snps = NULL,
                                                 header = TRUE,
                                                 snp_col = "markername",
                                                 beta_col = "beta",
                                                 se_col = "se_dgc",
                                                 eaf_col = "effect_allele_freq",
                                                 effect_allele_col = "effect_allele",
                                                 other_allele_col = "noneffect_allele",
                                                 pval_col = "p_dgc",
                                                 chr_col = "chr",
                                                 pos_col = "bp_hg19"))

colnames(cad) <- gsub(pattern = ".exposure","",colnames(cad))
cad$trait <- "coronary artery disease"
cad[,c("exposure","id","mr_keep","pval_origin")] <- NULL

data.table::fwrite(cad,"data/cad.txt", row.names = FALSE)
