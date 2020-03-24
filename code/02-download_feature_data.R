rm(list=setdiff(ls(), keep))
graphics.off()

# Download data ----------------------------------------------------------------

for (i in 1:nrow(features)) {
  download.file(features$source[i],paste0(path_features_raw,features$filename[i]),method = "wget", quiet = TRUE)
}

# Unzip files ------------------------------------------------------------------

for (i in 1:nrow(features)) {
  if (features$zip[i] %in% c(".bgz",".gz")) {
    system(paste0("gunzip -c ",path_features_raw,features$filename[i]," > ",path_features_processed,features$trait[i],features$ext[i]))
  } else if (features$zip[i] %in% c(".zip")) {
    system(paste0("unzip -c ",path_features_raw,features$filename[i]," > ",path_features_processed,features$trait[i],features$ext[i]))
  } else if (is.na(features$zip[i])) {
    system(paste0("cp ",path_features_raw,features$filename[i]," ",path_features_processed,features$trait[i],features$ext[i]))
  }
}