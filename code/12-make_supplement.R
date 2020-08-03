rm(list=ls())
graphics.off()

# Start excel workbook ----------------------------------------------------------

wb <- openxlsx::createWorkbook()

# Supplementary Table 1 - Features ----------------------------------

ST1 <- data.table::fread("raw/gwas.csv", data.table = FALSE)
openxlsx::addWorksheet(wb, "ST1")
openxlsx::writeData(wb, "ST1", ST1)

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