rm(list=ls())
graphics.off()

# Start excel workbook ----------------------------------------------------------

wb <- openxlsx::createWorkbook()

# Supplementary Table 1 - Features ----------------------------------

ST1 <- data.table::fread("raw/gwas.csv", data.table = FALSE)
openxlsx::addWorksheet(wb, "ST1")
openxlsx::writeData(wb, "ST1", ST1)

# Supplementary Table 2 - Instruments ----------------------------------

ST2 <- data.table::fread("data/instruments.csv", data.table = FALSE)
openxlsx::addWorksheet(wb, "ST2")
openxlsx::writeData(wb, "ST2", ST2)

# Save excel workbook ----------------------------------------------------------

openxlsx::saveWorkbook(wb, file = "output/SupplementaryTables.xlsx", overwrite = TRUE)

# Render markdown containing supplementary figures -----------------------------

rmarkdown::render('code/Rmd-supplementary_figures.Rmd', 
                  output_file =  "SupplementaryFigures.pdf", 
                  output_dir = 'output/')