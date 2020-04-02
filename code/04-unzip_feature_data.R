rm(list=setdiff(ls(), keep))
graphics.off()

# Load features data -----------------------------------------------------------

features <- data.table::fread("data/features_sources.csv",
                              stringsAsFactors = FALSE,
                              data.table = FALSE)

# Unzip files ------------------------------------------------------------------

files <- list.files(path = path_features_raw)
files <- c(files[41:43],files[45])

for (i in files) {
  if (grepl(".bgz",i)) {
    system(paste0("gunzip -c ",path_features_raw,i," > ",path_features_processed,features[features$filename==i,]$trait,features[features$filename==i,]$ext))
  } else if (grepl(".gz",i)) {
    system(paste0("gunzip -c ",path_features_raw,i," > ",path_features_processed,features[features$filename==i,]$trait,features[features$filename==i,]$ext))
  } else if (grepl(".zip",i)) {
    system(paste0("unzip -c ",path_features_raw,i," > ",path_features_processed,features[features$filename==i,]$trait,features[features$filename==i,]$ext))
  } else {
    system(paste0("cp ",path_features_raw,i," ",path_features_processed,features[features$filename==i,]$trait,features[features$filename==i,]$ext))
  }
}