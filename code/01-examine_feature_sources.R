rm(list=setdiff(ls(), keep))
graphics.off()

# Load feature data ------------------------------------------------------------

features <- readxl::read_xlsx("raw/feature_sources.xlsx",sheet = "All")

# Determine file name ----------------------------------------------------------

features$filename <-  ifelse(features$quick_download==TRUE,gsub(".*/","",features$source),NA)

# Determine file zip -----------------------------------------------------------

features$zip <-  ifelse(features$quick_download==TRUE,paste0(".",sub('.*\\.', '', features$filename)),NA)
features$zip <- ifelse(features$zip %in% c(".zip",".gz",".bgz"),features$zip,NA)

# Determine file extension -----------------------------------------------------

features$ext <- ifelse(features$quick_download==TRUE,paste0(".",gsub("\\..*","",gsub("^.*?\\.","",features$filename))),NA)
features$ext <- ifelse(features$consortium=="Lu",paste0(".",gsub(".*\\.","",gsub("^.*?\\.","",features$filename))),features$ext)
features$ext <- ifelse(features$consortium=="Neale",paste0(".",gsub(".*\\.","",gsub(".bgz","",features$filename[51]))),features$ext)

# Update adiponectin -----------------------------------------------------------

features[features$trait=="adiponectin",]$filename <- "adiponectin.txt"
features[features$trait=="adiponectin",]$zip <- NA
features[features$trait=="adiponectin",]$ext <- ".txt"

# Update liver enzymes ---------------------------------------------------------

enzymes <- features[features$consortium=="LOLIPOP",]$trait

for (i in enzymes) {
  features[features$trait==i,]$filename <- paste0(toupper(i),".txt.gz")
  features[features$trait==i,]$zip <- ".gz"
  features[features$trait==i,]$ext <- ".txt"
}


# Update CRP -----------------------------------------------------------

features[features$trait=="crp",]$filename <- "crp.txt"
features[features$trait=="crp",]$zip <- NA
features[features$trait=="crp",]$ext <- ".txt"

# Save -------------------------------------------------------------------------

data.table::fwrite(features,"data/features_sources.csv")