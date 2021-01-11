rm(list=ls())
graphics.off()

# Load results 

df <- data.table::fread("output/results_summary.csv", data.table = FALSE)

# Start excel workbook ----------------------------------------------------------

cad <- openxlsx::createWorkbook()
pad <- openxlsx::createWorkbook()

features <- setdiff(df$feature,c("t2d","pad","cad"))

for (i in features) {
  
  tmp <- df[df$feature==i,]
  
  tmp <- tidyr::pivot_longer(tmp, 
                             cols = colnames(tmp)[2:49], 
                             names_to = "name", 
                             values_to = "value")
  
  tmp$exposure_outcome <- gsub("\\..*","",tmp$name)
  tmp$stat <- gsub(".*\\.","",tmp$name)
  tmp$analysis <- substr(tmp$name,nchar(tmp$exposure_outcome)+2,nchar(tmp$name)-1-nchar(tmp$stat))
  tmp[,c("feature","name")] <- NULL
  tmp$analysis <- ifelse(tmp$analysis=="mvmr","mvmr - direct",tmp$analysis)
  tmp$analysis <- ifelse(tmp$analysis=="indir","mvmr - indirect",tmp$analysis)
  
  tmp <- tidyr::pivot_wider(tmp, names_from = "stat", values_from = "value")
  
  tmp_cad <- tmp[!grepl("pad",tmp$exposure_outcome),]
  openxlsx::addWorksheet(cad, i)
  openxlsx::writeData(cad, i, tmp_cad)
  
  tmp_pad <- tmp[!grepl("cad",tmp$exposure_outcome),]
  openxlsx::addWorksheet(pad, i)
  openxlsx::writeData(pad, i, tmp_pad)
  
}

openxlsx::saveWorkbook(cad, file = "output/results_summary_cad.xlsx", overwrite = TRUE)
openxlsx::saveWorkbook(pad, file = "output/results_summary_pad.xlsx", overwrite = TRUE)