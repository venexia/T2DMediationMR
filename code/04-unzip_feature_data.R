rm(list=setdiff(ls(), keep))
graphics.off()

# Load features data -----------------------------------------------------------

features <- data.table::fread("data/features_sources.csv",
                              stringsAsFactors = FALSE,
                              data.table = FALSE)

# Unzip files ------------------------------------------------------------------

for (i in 1:nrow(features)) {
  if (grepl(".bgz",features$filename[i])) {
    system(paste0("gunzip -c ",path_features_raw,features$filename[i]," > ",path_features_processed,features$trait[i],features$ext[i]))
  } else if (grepl(".gz",features$filename[i])) {
    system(paste0("gunzip -c ",path_features_raw,features$filename[i]," > ",path_features_processed,features$trait[i],features$ext[i]))
  } else if (grepl(".zip",features$filename[i])) {
    system(paste0("unzip -c ",path_features_raw,features$filename[i]," > ",path_features_processed,features$trait[i],features$ext[i]))
  } else {
    system(paste0("cp ",path_features_raw,features$filename[i]," ",path_features_processed,features$trait[i],features$ext[i]))
  }
}