# exp <- data.table::fread("raw/gwas-hba1c.txt", data.table = FALSE)
# 
# exp <- exp[exp$pval < 5e-8,]
# 
# exp <- TwoSampleMR::format_data(exp,
#                                 snp_col = "snp",
#                                 beta_col = "effect",
#                                 se_col = "stderr",
#                                 pval_col = "pvalue")
# 
# exp <- TwoSampleMR::clump_data(exp)
# 
# map <- data.frame(rbind(c("rs1046896","17:80685533"),
#              c("rs1387153","11:92673828"),
#              c("rs16926246","10:71093392"),
#              c("rs1799884","7:44229068"),
#              c("rs1800562","6:26093141"),
#              c("rs2779116","1:158585415"),
#              c("rs4737009","8:41630405"),
#              c("rs552976","2:169791438"),
#              c("rs6474359","8:41549194"),
#              c("rs7998202","13:113331868"),
#              c("rs855791","22:37462936")))
# 
# colnames(map) <- c("SNP","chr_pos")

exp <- data.table::fread("~/Downloads/HbA1c_METAL_European.txt.gz", data.table = FALSE)

exp <- exp[exp$pvalue < 5e-8,]

exp <- TwoSampleMR::format_data(exp,
                                type = "exposure",
                                beta_col = "beta",
                                chr_col = "chr",
                                pos_col = "pos",
                                effect_allele_col = "effect_allele",
                                other_allele_col = "other_allele",
                                eaf_col = "eaf_hapmap_CEU",
                                se_col = "stderr",
                                pval_col = "pvalue",
                                snp_col = "snp")

exp <- TwoSampleMR::clump_data(exp)

map <- unique(exp[,c("SNP","chr.exposure","pos.exposure")])
map$chr_pos <- paste(map$chr, map$pos, sep = ":")
map[,c("chr.exposure","pos.exposure")] <- NULL

out <- data.table::fread("raw/gwas-t2d.txt", data.table = FALSE)
out$chr_pos <- out$SNP
out$SNP <- NULL

out <- merge(out, map, by = "chr_pos")

out <- TwoSampleMR::format_data(out,
                                type = "outcome",
                                beta_col = "Beta",
                                chr_col = "Chr",
                                pos_col = "Pos",
                                effect_allele_col = "EA",
                                other_allele_col = "NEA",
                                eaf_col = "EAF",
                                se_col = "SE",
                                pval_col = "Pvalue",
                                snp_col = "SNP")

dat <- TwoSampleMR::harmonise_data(exp,out)

res <- TwoSampleMR::mr(dat)
