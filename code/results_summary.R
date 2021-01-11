rm(list=ls())
graphics.off()

# Load results 

uvmr <- data.table::fread("output/results.csv", data.table = FALSE)
uvmr <- uvmr[uvmr$method %in% c("Inverse variance weighted","Wald ratio"),]

mvmr <- data.table::fread("output/mvmr_results.csv", data.table = FALSE)
mvmr <- mvmr[mvmr$effect=="direct",]

# Create master dataframe

df <- data.frame(feature = unique(uvmr$exposure),
                 stringsAsFactors = FALSE)

# Add UVMR: Feature > T2D

tmp <- uvmr[uvmr$outcome=="t2d",c("exposure","b","se","pval","nsnp")]
colnames(tmp) <- c("feature","feature_t2d.uvmr.est","feature_t2d.uvmr.se","feature_t2d.uvmr.p","feature_t2d.uvmr.nsnp")
df <- merge(df, tmp, all.x = TRUE, by = "feature")

# Add UVMR: T2D > Feature

tmp <- uvmr[uvmr$exposure=="t2d",c("outcome","b","se","pval","nsnp")]
colnames(tmp) <- c("feature","t2d_feature.uvmr.est","t2d_feature.uvmr.se","t2d_feature.uvmr.p","t2d_feature.uvmr.nsnp")
df <- merge(df, tmp, all.x = TRUE, by = "feature")

# Add UVMR: Feature > CAD

tmp <- uvmr[uvmr$outcome=="cad",c("exposure","b","se","pval","nsnp")]
colnames(tmp) <- c("feature","feature_cad.uvmr.est","feature_cad.uvmr.se","feature_cad.uvmr.p","feature_cad.uvmr.nsnp")
df <- merge(df, tmp, all.x = TRUE, by = "feature")

# Add UVMR: Feature > PAD

tmp <- uvmr[uvmr$outcome=="pad",c("exposure","b","se","pval","nsnp")]
colnames(tmp) <- c("feature","feature_pad.uvmr.est","feature_pad.uvmr.se","feature_pad.uvmr.p","feature_pad.uvmr.nsnp")
df <- merge(df, tmp, all.x = TRUE, by = "feature")

# Add MVMR: Feature > CAD

tmp <- mvmr[mvmr$outcome=="cad" & mvmr$exposure!="t2d",c("exposure","estimate","se","pval")]
tmp$nsnp <- NA
colnames(tmp) <- c("feature","feature_cad.mvmr.est","feature_cad.mvmr.se","feature_cad.mvmr.p","feature_cad.mvmr.nsnp")
df <- merge(df, tmp, all.x = TRUE, by = "feature")

# Add MVMR: Feature > PAD

tmp <- mvmr[mvmr$outcome=="pad" & mvmr$exposure!="t2d",c("exposure","estimate","se","pval")]
tmp$nsnp <- NA
colnames(tmp) <- c("feature","feature_pad.mvmr.est","feature_pad.mvmr.se","feature_pad.mvmr.p","feature_pad.mvmr.nsnp")
df <- merge(df, tmp, all.x = TRUE, by = "feature")

# Add MVMR: T2D > CAD

tmp <- mvmr[mvmr$outcome=="cad" & mvmr$exposure=="t2d",c("analysis","estimate","se","pval")]
tmp$nsnp <- NA
tmp$analysis <- gsub("_cad","",tmp$analysis)
colnames(tmp) <- c("feature","t2d_cad.mvmr.est","t2d_cad.mvmr.se","t2d_cad.mvmr.p","t2d_cad.mvmr.nsnp")
df <- merge(df, tmp, all.x = TRUE, by = "feature")

# Add MVMR: T2D > PAD

tmp <- mvmr[mvmr$outcome=="pad" & mvmr$exposure=="t2d",c("analysis","estimate","se","pval")]
tmp$nsnp <- NA
tmp$analysis <- gsub("_pad","",tmp$analysis)
colnames(tmp) <- c("feature","t2d_pad.mvmr.est","t2d_pad.mvmr.se","t2d_pad.mvmr.p","t2d_pad.mvmr.nsnp")
df <- merge(df, tmp, all.x = TRUE, by = "feature")

# Add indirect: Feature > CAD

df$feature_cad.indir.est <- df$feature_t2d.uvmr.est*df$t2d_cad.mvmr.est
df$feature_cad.indir.se <- sqrt(df$feature_t2d.uvmr.se^2 + df$t2d_cad.mvmr.se^2)
df$feature_cad.indir.p <- NA  
df$feature_cad.indir.nsnp <- NA

# Add indirect: Feature > PAD

df$feature_pad.indir.est <- df$feature_t2d.uvmr.est*df$t2d_pad.mvmr.est
df$feature_pad.indir.se <- sqrt(df$feature_t2d.uvmr.se^2 + df$t2d_pad.mvmr.se^2)
df$feature_pad.indir.p <- NA  
df$feature_pad.indir.nsnp <- NA

# Add indirect: T2D > CAD

df$t2d_cad.indir.est <- df$t2d_feature.uvmr.est*df$feature_cad.mvmr.est
df$t2d_cad.indir.se <- sqrt(df$t2d_feature.uvmr.se^2 + df$feature_cad.mvmr.se^2)
df$t2d_cad.indir.p <- NA  
df$t2d_cad.indir.nsnp <- NA

# Add indirect: T2D > PAD

df$t2d_pad.indir.est <- df$t2d_feature.uvmr.est*df$feature_pad.mvmr.est
df$t2d_pad.indir.se <- sqrt(df$t2d_feature.uvmr.se^2 + df$feature_pad.mvmr.se^2)
df$t2d_pad.indir.p <- NA  
df$t2d_pad.indir.nsnp <- NA

# Save

data.table::fwrite(df,"output/results_summary.csv")