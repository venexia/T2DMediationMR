rm(list=setdiff(ls(), keep))
graphics.off()

# Load features data -----------------------------------------------------------

features <- data.table::fread("data/features_sources.csv",
                              colClasses = c(rep("character",6),"numeric",rep("character",19)),
                              stringsAsFactors = FALSE,
                              data.table = FALSE)

features <- features[features$quick_download==TRUE,]

# Download data ----------------------------------------------------------------

for (i in 1:nrow(features)) {
  download.file(features$source[i],paste0(path_features_raw,features$filename[i]),method = "wget", quiet = TRUE)
}