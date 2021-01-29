rm(list=ls())
graphics.off()

# Start excel workbook ----------------------------------------------------------

wb <- openxlsx::createWorkbook()

# Supplementary Table 1 - Features ---------------------------------------------

ST1 <- data.table::fread("raw/gwas.csv", data.table = FALSE)
ST1 <- ST1[ST1$source!="exclude_feature",]
ST1$feature <- ST1$source!="outcome"
ST1$source <- NULL
openxlsx::addWorksheet(wb, "ST1")
openxlsx::writeData(wb, "ST1", ST1)

# Supplementary Table 2 - Univariate estimates using alternative methods -------

ST2 <- data.table::fread("output/results.csv", data.table = FALSE)

ST2 <- merge(ST2, ST1[,c("trait","trait_long")], by.x = "exposure", by.y = "trait")
ST2$exposure <- ST2$trait_long
ST2$trait_long <- NULL
ST2 <- merge(ST2, ST1[,c("trait","trait_long")], by.x = "outcome", by.y = "trait")
ST2$outcome <- ST2$trait_long
ST2$trait_long <- NULL

ST2 <- ST2[,c("exposure","outcome","method","nsnp","b","se","pval")]
colnames(ST2) <- c("exposure","outcome","method","nsnp","estimate","se","pvalue")

ST2$method <- ifelse(ST2$method %in% c("Wald ratio","Inverse variance weighted"),"Main",ST2$method)
ST2$method <- gsub(" ","_",ST2$method)

ST2 <- tidyr::pivot_wider(ST2, 
                          names_from = "method",
                          names_glue = "{method}_{.value}",
                          values_from = c("estimate","se","pvalue"),
                          values_fill = NA)

ST2 <- ST2[!(ST2$exposure %in% c("coronary heart disease","peripheral artery disease")),]
ST2 <- ST2[,c("exposure","outcome","nsnp",
              paste0("Main",c("_estimate","_se","_pvalue")),
              paste0("MR_Egger",c("_estimate","_se","_pvalue")),
              paste0("Simple_mode",c("_estimate","_se","_pvalue")),
              paste0("Weighted_mode",c("_estimate","_se","_pvalue")),
              paste0("Weighted_median",c("_estimate","_se","_pvalue")))]

ST2 <- ST2[order(ST2$exposure,ST2$outcome),]

openxlsx::addWorksheet(wb, "ST2")
openxlsx::writeData(wb, "ST2", ST2)

# Supplementary Table 3 - two step MR estimates --------------------------------

ST3 <- data.table::fread("output/twostep_results.csv", data.table = FALSE)
ST3 <- ST3[ST3$effect %in% c("direct","indirect"),c("analysis","effect","exposure","mediator","outcome","estimate","se","Qstat","Qpval","Qdf","condF")]
ST3 <- merge(ST3, ST1[,c("trait","trait_long")], by.x = "exposure", by.y = "trait")
ST3$exposure <- ST3$trait_long
ST3$trait_long <- NULL
ST3 <- merge(ST3, ST1[,c("trait","trait_long")], by.x = "mediator", by.y = "trait")
ST3$mediator <- ST3$trait_long
ST3$trait_long <- NULL
ST3 <- merge(ST3, ST1[,c("trait","trait_long")], by.x = "outcome", by.y = "trait")
ST3$outcome <- ST3$trait_long
ST3$trait_long <- NULL
ST3 <- ST3[,c("exposure","mediator","outcome","effect","estimate","se","Qstat","Qpval","Qdf","condF")]
openxlsx::addWorksheet(wb, "ST3")
openxlsx::writeData(wb, "ST3", ST3)

# Supplementary Table 4 - Egger intercept test ---------------------------------

ST4 <- data.table::fread("output/plei.csv", data.table = FALSE)
ST4 <- merge(ST4, ST1[,c("trait","trait_long")], by.x = "exposure", by.y = "trait", all.y = TRUE)
ST4$exposure <- ST4$trait_long
ST4$trait_long <- NULL
ST4 <- merge(ST4, ST1[,c("trait","trait_long")], by.x = "outcome", by.y = "trait")
ST4$outcome <- ST4$trait_long
ST4$trait_long <- NULL
ST4 <- ST4[,c("exposure","outcome","egger_intercept","se","pval")]
ST4 <- na.omit(ST4)
ST4 <- ST4[!(ST4$exposure %in% c("coronary heart disease","peripheral artery disease")),]
openxlsx::addWorksheet(wb, "ST4")
openxlsx::writeData(wb, "ST4", ST4)

# Supplementary Table 5 - I squared statistics for MVMR ------------------------

ST5 <- data.table::fread("output/results.csv", data.table = FALSE)
ST5 <- ST5[,c("exposure","outcome","nsnp","Isq")]
ST5 <- unique(ST5)
ST5 <- merge(ST5, ST1[,c("trait","trait_long")], by.x = "exposure", by.y = "trait", all.x = TRUE)
ST5$exposure <- ST5$trait_long
ST5$trait_long <- NULL
ST5 <- merge(ST5, ST1[,c("trait","trait_long")], by.x = "outcome", by.y = "trait", all.x = TRUE)
ST5$outcome <- ST5$trait_long
ST5$trait_long <- NULL
ST5 <- ST5[,c("exposure","outcome","Isq")]
ST5 <- na.omit(ST5)
ST5 <- ST5[!(ST5$exposure %in% c("coronary heart disease","peripheral artery disease")),]
openxlsx::addWorksheet(wb, "ST5")
openxlsx::writeData(wb, "ST5", ST5)

# Save excel workbook ----------------------------------------------------------

openxlsx::saveWorkbook(wb, file = "output/SupplementaryTables.xlsx", overwrite = TRUE)

# Supplementary Content - Instruments ----------------------------------

wb2 <- openxlsx::createWorkbook()
instruments <- data.table::fread("data/instruments_all.txt", data.table = FALSE)
for (i in unique(instruments$exposure)) {
  openxlsx::addWorksheet(wb2, i)
  openxlsx::writeData(wb2, i, instruments[instruments$exposure==i,])
}
openxlsx::saveWorkbook(wb2, file = "output/Instruments.xlsx", overwrite = TRUE)

# Render markdown containing supplementary figures -----------------------------

rmarkdown::render('code/Rmd-supplementary_figures.Rmd', 
                  output_file =  "SupplementaryFigures.pdf", 
                  output_dir = 'output/')