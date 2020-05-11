rm(list=setdiff(ls(), keep))
graphics.off()

# Clean CAD GWAS ---------------------------------------------------------------

cad <- data.table::fread("raw/cad.add.160614.website.txt",
                         stringsAsFactors = FALSE,
                         data.table = FALSE)

cad$trait <- "coronary artery disease"

cad <- cad[,c("markername","chr","bp_hg19","effect_allele","noneffect_allele","effect_allele_freq","beta","se_dgc","p_dgc","trait")]

colnames(cad) <- c("SNP","chr","pos","effect_allele","other_allele","eaf","beta","se","pval","trait")

data.table::fwrite(cad,"data/cad.txt", row.names = FALSE)

# Clean PAD GWAS ---------------------------------------------------------------

pad <- data.table::fread(paste0(path_pad_gwas,"CLEANED.MVP.EUR.PAD.results.anno.csv"),
                         stringsAsFactors = FALSE,
                         data.table = FALSE)

pad$trait <- "peripheral artery disease"

pad <- pad[,c("ID","CHROM","POS","EFFECT_ALLELE","OTHER_ALLELE","EAF","BETA","SE","PVAL","trait")]

colnames(pad) <- c("SNP","chr","pos","effect_allele","other_allele","eaf","beta","se","pval","trait")

data.table::fwrite(pad,"data/pad.txt", row.names = FALSE)