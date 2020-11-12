
gwas_logistic <- data.table::fread("data/gwas-pad.txt", data.table = FALSE)
gwas_linear <- data.table::fread("raw/PAD.glm.linear", data.table = FALSE)

gwas_logistic$ID1 <- paste(gwas_logistic$chr, gwas_logistic$pos, 
                          gwas_logistic$effect_allele, gwas_logistic$other_allele, sep = ":")

gwas_logistic$ID2 <- paste(gwas_logistic$chr, gwas_logistic$pos, 
                           gwas_logistic$other_allele, gwas_logistic$effect_allele, sep = ":")

gwas_logistic2 <- gwas_logistic[(gwas_logistic$ID1 %in% gwas_linear$ID | gwas_logistic$ID2 %in% gwas_linear$ID),]

tmp_gwas_linear <- gwas_linear[,c("CHR","POS","A1","BETA")]
colnames(tmp_gwas_linear) <- c("CHR","POS","A1_linear","beta_linear")

tmp_gwas_logistic <- gwas_logistic2[,c("chr","pos","effect_allele","beta")]
colnames(tmp_gwas_logistic) <- c("CHR","POS","A1_logistic","beta_logistic")

df <- merge(tmp_gwas_linear, tmp_gwas_logistic, by = c("CHR","POS"))