rm(list=ls())
graphics.off()

df <-  data.table::fread("output/mvmr_results_indirect.csv",
                         stringsAsFactors = FALSE,
                         data.table = FALSE)

tmp <- tidyr::pivot_wider(df, names_from = "effect", values_from = c("estimate","se"))

tmp$estimate_total_2stepMVMR <- tmp$estimate_indirect_2stepMVMR + tmp$estimate_direct
tmp$se_total_2stepMVMR <- NA

tmp$estimate_total_2stepUVMR <- tmp$estimate_indirect_2stepUVMR + tmp$estimate_direct
tmp$se_total_2stepUVMR <- NA

tmp <- tmp[,c("analysis","exposure","outcome","estimate_total_2stepMVMR","se_total_2stepMVMR","estimate_total_2stepUVMR","se_total_2stepUVMR")]
tmp <- tmp[!is.na(tmp$estimate_total_2stepMVMR),]
tmp <- tidyr::pivot_longer(tmp, cols = setdiff(colnames(tmp),c("analysis","exposure","outcome")))
tmp$estimate_se <- gsub("_.*","",tmp$name)
tmp$effect <- gsub("estimate_","",tmp$name)
tmp$effect <- gsub("se_","",tmp$effect)
tmp$name <- NULL
tmp <- tidyr::pivot_wider(tmp, names_from = "estimate_se", values_from = "value")

df <- rbind(df,tmp)