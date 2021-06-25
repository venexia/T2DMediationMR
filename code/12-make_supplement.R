rm(list=ls())
graphics.off()

# Start excel workbook ----------------------------------------------------------

wb <- openxlsx::createWorkbook()

# Supplementary Table 1 - Features ---------------------------------------------

risk_factors <- data.table::fread("data/risk_factors.csv", data.table = FALSE)
outcomes <- data.table::fread("raw/outcomes.csv", data.table = FALSE)
ST1 <- rbind(risk_factors,outcomes)

openxlsx::addWorksheet(wb, "ST1")
openxlsx::writeData(wb, "ST1", ST1)

# Supplementary Table 2 - Univariate estimates using alternative methods -------

fwd <- data.table::fread("output/results_fwd.csv", data.table = FALSE)
bkwd <- data.table::fread("output/results_bkwd.csv", data.table = FALSE)
ST2 <- rbind(fwd,bkwd)
ST2 <- ST2[ST2$exposure!=ST2$outcome,]

ST2 <- merge(ST2, ST1[,c("id","trait")], by.x = "exposure", by.y = "id", all.x = TRUE)
ST2 <- dplyr::rename(ST2, exposure_id = exposure, exposure = trait)

ST2 <- merge(ST2, ST1[,c("id","trait")], by.x = "outcome", by.y = "id")
ST2 <- dplyr::rename(ST2, outcome_id = outcome, outcome = trait)

ST2 <- ST2[,c("exposure","exposure_id","outcome","outcome_id","method","nsnp","b","se","pval")]
colnames(ST2) <- c("exposure","exposure_id","outcome","outcome_id","method","nsnp","estimate","se","pvalue")

ST2$method <- ifelse(ST2$method %in% c("Wald ratio","Inverse variance weighted"),"Main",ST2$method)
ST2$method <- gsub(" ","_",ST2$method)

ST2 <- tidyr::pivot_wider(ST2, 
                          names_from = "method",
                          names_glue = "{method}_{.value}",
                          values_from = c("estimate","se","pvalue"),
                          values_fill = NA)

ST2 <- ST2[!(ST2$exposure %in% c("Coronary heart disease","Peripheral artery disease")),]
ST2 <- ST2[,c("exposure","outcome","nsnp",
              paste0("Main",c("_estimate","_se","_pvalue")),
              paste0("MR_Egger",c("_estimate","_se","_pvalue")),
              paste0("Simple_mode",c("_estimate","_se","_pvalue")),
              paste0("Weighted_mode",c("_estimate","_se","_pvalue")),
              paste0("Weighted_median",c("_estimate","_se","_pvalue")))]

ST2 <- ST2[order(ST2$exposure,ST2$outcome),]

openxlsx::addWorksheet(wb, "ST2")
openxlsx::writeData(wb, "ST2", ST2)

# Supplementary Table 3 - evidence summary -------------------------------------

ST3 <- data.table::fread("output/evidence_summary.csv")
colnames(ST3) <- gsub("_adjust","_fdr",colnames(ST3))
ST3 <- dplyr::rename(ST3, 
              rf = trait,
              rf_id = id,
              rf_t2d_pval = rf_t2d,
              t2d_rf_pval = t2d_rf,
              rf_pad_pval = rf_pad,
              rf_cad_pval = rf_cad)

openxlsx::addWorksheet(wb, "ST3")
openxlsx::writeData(wb, "ST3", ST3)

# Supplementary Table 4 - two step MR estimates --------------------------------

ST4 <- data.table::fread("output/twostep_results.csv", data.table = FALSE)
ST4 <- ST4[ST4$effect %in% c("direct","indirect"),c("effect","exposure","mediator","outcome","estimate","se","Qstat","Qpval","Qdf","condF")]

ST4 <- merge(ST4, ST1[,c("id","trait")], by.x = "exposure", by.y = "id")
ST4 <- dplyr::rename(ST4, exposure_id = exposure, exposure = trait)

ST4 <- merge(ST4, ST1[,c("id","trait")], by.x = "mediator", by.y = "id")
ST4 <- dplyr::rename(ST4, mediator_id = mediator, mediator = trait)

ST4 <- merge(ST4, ST1[,c("id","trait")], by.x = "outcome", by.y = "id")
ST4 <- dplyr::rename(ST4, outcome_id = outcome, outcome = trait)

ST4 <- ST4[,c("exposure","exposure_id","mediator","mediator_id","outcome","outcome_id","effect","estimate","se","Qstat","Qpval","Qdf","condF")]

ST4 <- ST4[order(ST4$exposure,ST4$mediator,ST4$outcome,ST4$effect),]
openxlsx::addWorksheet(wb, "ST4")
openxlsx::writeData(wb, "ST4", ST4)

# Supplementary Table 5 - Egger intercept test ---------------------------------

plei_fwd <- data.table::fread("output/plei_fwd.csv", data.table = FALSE)
plei_bkwd <- data.table::fread("output/plei_bkwd.csv", data.table = FALSE)
ST5 <- rbind(plei_fwd, plei_bkwd)

ST5 <- merge(ST5, ST1[,c("id","trait")], by.x = "exposure", by.y = "id", all.x = TRUE)
ST5 <- dplyr::rename(ST5, exposure_id = exposure, exposure = trait)

ST5 <- merge(ST5, ST1[,c("id","trait")], by.x = "outcome", by.y = "id")
ST5 <- dplyr::rename(ST5, outcome_id = outcome, outcome = trait)

ST5 <- ST5[,c("exposure","exposure_id","outcome","outcome_id","egger_intercept","se","pval")]
ST5 <- na.omit(ST5)
ST5 <- ST5[!(ST5$exposure %in% c("coronary heart disease","peripheral artery disease")),]

ST5 <- ST5[order(ST5$exposure,ST5$outcome),]
openxlsx::addWorksheet(wb, "ST5")
openxlsx::writeData(wb, "ST5", ST5)

# Supplementary Table 6 - I squared statistics for MVMR ------------------------

fwd <- data.table::fread("output/results_fwd.csv", 
                         select = c("exposure","outcome","nsnp","Isq"),
                         data.table = FALSE)

bkwd <- data.table::fread("output/results_bkwd.csv", 
                          select = c("exposure","outcome","nsnp","Isq"),
                          data.table = FALSE)

ST6 <- rbind(fwd,bkwd)
ST6 <- ST6[ST6$exposure!=ST6$outcome,]
ST6 <- unique(ST6)

ST6 <- merge(ST6, ST1[,c("id","trait")], by.x = "exposure", by.y = "id", all.x = TRUE)
ST6 <- dplyr::rename(ST6, exposure_id = exposure, exposure = trait)

ST6 <- merge(ST6, ST1[,c("id","trait")], by.x = "outcome", by.y = "id")
ST6 <- dplyr::rename(ST6, outcome_id = outcome, outcome = trait)

ST6 <- ST6[,c("exposure","exposure_id","outcome","outcome_id","Isq")]
ST6 <- na.omit(ST6)
ST6 <- ST6[!(ST6$exposure %in% c("coronary heart disease","peripheral artery disease")),]

ST6 <- ST6[order(ST6$exposure,ST6$outcome),]
openxlsx::addWorksheet(wb, "ST6")
openxlsx::writeData(wb, "ST6", ST6)

# Save excel workbook ----------------------------------------------------------

openxlsx::saveWorkbook(wb, file = "output/SupplementaryTables.xlsx", overwrite = TRUE)

# Supplementary Content - Instruments ----------------------------------

wb2 <- openxlsx::createWorkbook()
instruments <- data.table::fread("data/instruments.csv", data.table = FALSE)
for (i in unique(instruments$exposure)) {
  openxlsx::addWorksheet(wb2, i)
  openxlsx::writeData(wb2, i, instruments[instruments$exposure==i,])
}
openxlsx::saveWorkbook(wb2, file = "output/Instruments.xlsx", overwrite = TRUE)

# Render markdown containing supplementary figures -----------------------------

rmarkdown::render('code/Rmd-supplementary_figures.Rmd', 
                  output_file =  "SupplementaryFigures.pdf", 
                  output_dir = 'output/')