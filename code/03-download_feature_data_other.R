rm(list=setdiff(ls(), keep))
graphics.off()

# # Adiponectin ------------------------------------------------------------------
# 
# adiponectin_files <- list.files(path = "raw/", 
#                                 pattern = "adipogen.discovery.eur_.meta_.public.release.part")
# 
# adiponectin <- NULL
# 
# for (file in adiponectin_files) {
#   tmp <- data.table::fread(paste0("raw/",file), data.table = FALSE, stringsAsFactors = FALSE)
#   adiponectin <- rbind(adiponectin,tmp)
# }
# 
# data.table::fwrite(adiponectin,paste0(path_features_raw,"adiponectin.txt"))

# Liver enzymes ----------------------------------------------------------------

features <- data.table::fread("data/features_sources.csv", data.table = FALSE, stringsAsFactors = FALSE)

enzymes <- features[features$consortium=="LOLIPOP",]$trait

for (i in enzymes) {
  system(paste0("cp raw/",toupper(i),".txt.gz ",path_features_raw,toupper(i),".txt.gz"))
}

# CRP ----------------------------------------------------------------

system(paste0("cp raw/GWAS_CRP_eu_dv1.CRP.1.txt ",path_features_raw,"crp.txt"))
