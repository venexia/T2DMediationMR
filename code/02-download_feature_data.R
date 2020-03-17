features <- readxl::read_xlsx("raw/feature_sources.xlsx",sheet = "All")
features <- features[features$quick_download==TRUE,]

# Determine file name ----------------------------------------------------------

features$filename <- gsub(".*/","",features$source)

# Determine file type ----------------------------------------------------------

# features$filetype <- ifelse(features$consortium=="Neale",".tsv.bgz",gsub("^[^.]*","",features$filename))
# 
# features$sv <- NA
# features$sv <- ifelse(features$filetype %in% c(".TBL.txt",".txt.gz",".txt"),".txt",features$sv)
# features$sv <- ifelse(features$filetype %in% c(".csv.zip"),".csv",features$sv)
# features$sv <- ifelse(features$filetype %in% c(".tsv.bgz"),".tsv",features$sv)
# 
# features$zip <- NA
# features$zip <- ifelse(grepl(".gz",features$filetype),".gz",features$zip)
# features$zip <- ifelse(grepl(".zip",features$filetype),".zip",features$zip)
# features$zip <- ifelse(grepl(".bgz",features$filetype),".bgz",features$zip)

# Download data ----------------------------------------------------------------

for (i in 1:nrow(features)) {
  
  system(paste0("curl -o ",path_features,features$filename[i]," ",features$source[i]))
  
  # if (!is.na(features$zip[i])) {
  #   system(paste0("unzip ",path_features,features$filename[i]))
  #   system(paste0("cp ",path_features,gsub(features$zip[i],"",features$filename[i])," ",path_features,features$trait[i],features$sv[i]))
  #   system(paste0("rm ",path_features,features$filename[i]))
  #   system(paste0("rm ",path_features,gsub(features$zip[i],"",features$filename[i])))
  # } else {
  #   system(paste0("cp ",path_features,features$filename[i]," ",path_features,features$trait[i],features$sv[i]))
  #   system(paste0("rm ",path_features,features$filename[i]))
  # }
  # 
  # if (features$sv[i]!=".txt") {
  #   tmp <- data.table::fread(paste0(path_features,features$trait[i],features$sv[i]))
  #   data.table::fwrite(tmp,paste0(path_features,features$trait[i],".txt"))
  #   system(paste0("rm ",path_features,features$trait[i],features$sv[i]))
  # }  
  
}