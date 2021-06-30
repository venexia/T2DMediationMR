rm(list=ls())
graphics.off()

# Start excel workbook ----------------------------------------------------------

wb <- openxlsx::createWorkbook()

c <- 0

# Supplementary Table - Risk factor---------------------------------------------

c <- c + 1

risk_factors <- data.table::fread("data/risk_factors_10snps.csv", data.table = FALSE)
outcomes <- data.table::fread("raw/outcomes.csv", data.table = FALSE)
gwas <- rbind(risk_factors,outcomes)
gwas$trait <- gsub("  automated.*","",gwas$trait)
gwas$trait <- ifelse(gwas$trait=="Alcohol intake frequency.","Alcohol intake frequency",gwas$trait)

openxlsx::addWorksheet(wb, paste0("ST",c))
openxlsx::writeData(wb, paste0("ST",c), gwas)

# Supplementary Table - Univariate estimates using alternative methods ---------

c <- c + 1

fwd <- data.table::fread("output/results_fwd.csv", data.table = FALSE)
bkwd <- data.table::fread("output/results_bkwd.csv", data.table = FALSE)
df <- rbind(fwd,bkwd)
df <- df[df$exposure!=df$outcome,]

df <- merge(df, gwas[,c("id","trait")], by.x = "exposure", by.y = "id", all.x = TRUE)
df <- dplyr::rename(df, exposure_id = exposure, exposure = trait)

df <- merge(df, gwas[,c("id","trait")], by.x = "outcome", by.y = "id")
df <- dplyr::rename(df, outcome_id = outcome, outcome = trait)

df <- df[,c("exposure","exposure_id","outcome","outcome_id","method","nsnp","b","se","pval")]
colnames(df) <- c("exposure","exposure_id","outcome","outcome_id","method","nsnp","estimate","se","pvalue")

df$method <- ifelse(df$method %in% c("Wald ratio","Inverse variance weighted"),"Main",df$method)
df$method <- gsub(" ","_",df$method)

df <- tidyr::pivot_wider(df, 
                          names_from = "method",
                          names_glue = "{method}_{.value}",
                          values_from = c("estimate","se","pvalue"),
                          values_fill = NA)

df <- df[!(df$exposure %in% c("Coronary heart disease","Peripheral artery disease")),]
df <- df[,c("exposure","outcome","nsnp",
              paste0("Main",c("_estimate","_se","_pvalue")),
              paste0("MR_Egger",c("_estimate","_se","_pvalue")),
              paste0("Simple_mode",c("_estimate","_se","_pvalue")),
              paste0("Weighted_mode",c("_estimate","_se","_pvalue")),
              paste0("Weighted_median",c("_estimate","_se","_pvalue")))]

df <- df[order(df$exposure,df$outcome),]

openxlsx::addWorksheet(wb, paste0("ST",c))
openxlsx::writeData(wb, paste0("ST",c), df)

# Supplementary Table - Egger intercept test -----------------------------------

c <- c + 1

plei_fwd <- data.table::fread("output/plei_fwd.csv", data.table = FALSE)
plei_bkwd <- data.table::fread("output/plei_bkwd.csv", data.table = FALSE)
df <- rbind(plei_fwd, plei_bkwd)

df <- merge(df, gwas[,c("id","trait")], by.x = "exposure", by.y = "id", all.x = TRUE)
df <- dplyr::rename(df, exposure_id = exposure, exposure = trait)

df <- merge(df, gwas[,c("id","trait")], by.x = "outcome", by.y = "id")
df <- dplyr::rename(df, outcome_id = outcome, outcome = trait)

df <- df[,c("exposure","exposure_id","outcome","outcome_id","egger_intercept","se","pval")]
df <- na.omit(df)
df <- df[!(df$exposure %in% c("coronary heart disease","peripheral artery disease")),]

df <- df[order(df$exposure,df$outcome),]
openxlsx::addWorksheet(wb, paste0("ST",c))
openxlsx::writeData(wb, paste0("ST",c), df)

# Supplementary Table - I squared statistics for MR Egger ----------------------

c <- c + 1

fwd <- data.table::fread("output/results_fwd.csv", 
                         select = c("exposure","outcome","nsnp","Isq"),
                         data.table = FALSE)

bkwd <- data.table::fread("output/results_bkwd.csv", 
                          select = c("exposure","outcome","nsnp","Isq"),
                          data.table = FALSE)

df <- rbind(fwd,bkwd)
df <- df[df$exposure!=df$outcome,]
df <- unique(df)

df <- merge(df, gwas[,c("id","trait")], by.x = "exposure", by.y = "id", all.x = TRUE)
df <- dplyr::rename(df, exposure_id = exposure, exposure = trait)

df <- merge(df, gwas[,c("id","trait")], by.x = "outcome", by.y = "id")
df <- dplyr::rename(df, outcome_id = outcome, outcome = trait)

df <- df[,c("exposure","exposure_id","outcome","outcome_id","Isq")]
df <- na.omit(df)
df <- df[!(df$exposure %in% c("coronary heart disease","peripheral artery disease")),]

df <- df[order(df$exposure,df$outcome),]
openxlsx::addWorksheet(wb, paste0("ST",c))
openxlsx::writeData(wb, paste0("ST",c), df)

# Supplementary Table - evidence summary ---------------------------------------

c <- c + 1

df <- data.table::fread("output/evidence_summary.csv")
colnames(df) <- gsub("_adjust","_fdr",colnames(df))
df <- dplyr::rename(df, 
                    rf = trait,
                    rf_id = id,
                    rf_t2d_pval = rf_t2d,
                    t2d_rf_pval = t2d_rf,
                    rf_pad_pval = rf_pad,
                    rf_cad_pval = rf_cad)

openxlsx::addWorksheet(wb, paste0("ST",c))
openxlsx::writeData(wb, paste0("ST",c), df)

# Supplementary Table - two step MR estimates ----------------------------------

c <- c + 1

df <- data.table::fread("output/twostep_results.csv", data.table = FALSE)
df <- df[df$effect %in% c("direct","indirect"),c("effect","exposure","mediator","outcome","estimate","se","Qstat","Qpval","Qdf","condF")]

df <- merge(df, gwas[,c("id","trait")], by.x = "exposure", by.y = "id")
df <- dplyr::rename(df, exposure_id = exposure, exposure = trait)

df <- merge(df, gwas[,c("id","trait")], by.x = "mediator", by.y = "id")
df <- dplyr::rename(df, mediator_id = mediator, mediator = trait)

df <- merge(df, gwas[,c("id","trait")], by.x = "outcome", by.y = "id")
df <- dplyr::rename(df, outcome_id = outcome, outcome = trait)

df <- df[,c("exposure","exposure_id","mediator","mediator_id","outcome","outcome_id","effect","estimate","se","Qstat","Qpval","Qdf","condF")]

df <- df[order(df$exposure,df$mediator,df$outcome,df$effect),]
openxlsx::addWorksheet(wb, paste0("ST",c))
openxlsx::writeData(wb, paste0("ST",c), df)

# Save excel workbook ----------------------------------------------------------

openxlsx::saveWorkbook(wb, file = "output/SupplementaryTables.xlsx", overwrite = TRUE)

# Supplementary Content - Instruments ----------------------------------

wb2 <- openxlsx::createWorkbook()
instruments <- data.table::fread("data/instruments.csv", data.table = FALSE)
instruments <- merge(instruments, gwas[,c("id","trait")], by.x = "exposure", by.y = "id")
instruments <- dplyr::rename(instruments, exposure_id = exposure, exposure = trait)

for (i in unique(instruments$exposure_id)) {
  openxlsx::addWorksheet(wb2, i)
  openxlsx::writeData(wb2, i, instruments[instruments$exposure_id==i,])
}
openxlsx::saveWorkbook(wb2, file = "output/Instruments.xlsx", overwrite = TRUE)

# Render markdown containing supplementary figures -----------------------------

rmarkdown::render('code/Rmd-supplementary_figures.Rmd', 
                  output_file =  "SupplementaryFigures.pdf", 
                  output_dir = 'output/')