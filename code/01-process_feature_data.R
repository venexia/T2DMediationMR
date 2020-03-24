rm(list=setdiff(ls(), keep))
graphics.off()

# Load feature data ------------------------------------------------------------

features <- readxl::read_xlsx("raw/feature_sources.xlsx",sheet = "All")
features <- features[features$quick_download==TRUE,]

# Determine file name ----------------------------------------------------------

features$filename <- gsub(".*/","",features$source)

# Decompose file name ----------------------------------------------------------

features$zip <- paste0(".",sub('.*\\.', '', features$filename))
features$zip <- ifelse(features$zip %in% c(".zip",".gz",".bgz"),features$zip,NA)

features$ext <- paste0(".",gsub("\\..*","",gsub("^.*?\\.","",features$filename)))
features$ext <- ifelse(features$consortium=="Lu",paste0(".",gsub(".*\\.","",gsub("^.*?\\.","",features$filename))),features$ext)
features$ext <- ifelse(features$consortium=="Neale",paste0(".",gsub(".*\\.","",gsub(".bgz","",features$filename[51]))),features$ext)

# Save enhanced features speadsheet --------------------------------------------

data.table::fwrite(features,"data/features_sources.csv")