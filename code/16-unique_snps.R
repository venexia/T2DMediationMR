df <-  data.table::fread("output/mvmr_results.csv",
                         stringsAsFactors = FALSE,
                         #select = c("analysis","exposure","outcome","effect","estimate","se"),
                         data.table = FALSE)

snps <- paste0(df$snps,collapse = ";")
snps <- strsplit(snps,";")
snps <- as.vector(snps[[1]])
snps <- unique(snps)