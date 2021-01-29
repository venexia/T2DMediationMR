rm(list=ls())
graphics.off()

# Set FDR threshold for features to progress to next analysis ------------------

threshold <- 0.05
outcomes <- c("t2d","pad","cad")

# Load univariate results ------------------------------------------------------

df <- data.table::fread("output/results.csv",
                        select = c("exposure", "outcome", "method","pval"),
                        data.table = FALSE)

df <- df[df$method %in% c("Inverse variance weighted","Wald ratio"),]

df$analysis <- ""

df$analysis <- ifelse(df$exposure=="t2d" & !(df$outcome %in% outcomes),
                      "t2d_feature",
                      df$analysis)

df$analysis <- ifelse(!(df$exposure %in% outcomes) & df$outcome=="t2d",
                      "feature_t2d",
                      df$analysis)

df$analysis <- ifelse(!(df$exposure %in% outcomes) & df$outcome=="cad",
                      "feature_cad",
                      df$analysis)

df$analysis <- ifelse(!(df$exposure %in% outcomes) & df$outcome=="pad",
                      "feature_pad",
                      df$analysis)

df <- df[df$analysis!="",]

df$feature <- ifelse(df$analysis=="t2d_feature",df$outcome,df$exposure)

df <- df[,c("feature","analysis","pval")]

df <- tidyr::pivot_wider(df, names_from = "analysis", values_from = "pval")

for (i in c("feature_t2d","t2d_feature","feature_pad","feature_cad")) {
  cols <- colnames(df)
  tmp <- data.frame(cols = df[,i])
  colnames(tmp) <- c("pval")
  df$BH <- stats::p.adjust(tmp$pval, method = "BH")
  df$evidence <- df$BH<threshold
  colnames(df) <- c(cols,paste0(i,c("_adjust","_evidence")))
}

df <- df[,c("feature","feature_t2d_evidence","t2d_feature_evidence","feature_pad_evidence","feature_cad_evidence")]
df <- na.omit(df)

data.table::fwrite(df,"output/evidence_summary.csv")