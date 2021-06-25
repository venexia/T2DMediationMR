rm(list=ls())
graphics.off()

# Set FDR threshold for risk factors to progress to next analysis ------------------

threshold <- 0.05

# Load master list -------------------------------------------------------------

master <- data.table::fread("data/risk_factors.csv", 
                            select = c("trait","id"),
                            data.table = FALSE)

# Load univariate results ------------------------------------------------------

df_fwd <- data.table::fread("output/results_fwd.csv",
                        select = c("exposure", "outcome", "method","pval"),
                        data.table = FALSE)

df_bkwd <- data.table::fread("output/results_bkwd.csv",
                            select = c("exposure", "outcome", "method","pval"),
                            data.table = FALSE)

df <- rbind(df_fwd, df_bkwd)

df <- df[df$outcome!="t2d_linear" & df$exposure!="t2d_linear",]
df <- df[df$method %in% c("Inverse variance weighted","Wald ratio"),]
df <- df[!(df$exposure %in% c("pad","cad")),]
df <- df[!(df$exposure=="t2d" & df$outcome %in% c("pad","cad","t2d")),]

df$analysis <- ""
df$analysis <- ifelse(df$exposure=="t2d","t2d_rf",df$analysis)
df$analysis <- ifelse(df$outcome=="t2d","rf_t2d",df$analysis)
df$analysis <- ifelse(df$outcome=="cad","rf_cad",df$analysis)
df$analysis <- ifelse(df$outcome=="pad","rf_pad",df$analysis)

df$id <- ifelse(df$analysis=="t2d_rf",df$outcome,df$exposure)
df <- df[,c("id","analysis","pval")]

df <- tidyr::pivot_wider(df, names_from = "analysis", values_from = "pval")

for (i in c("rf_t2d","t2d_rf","rf_pad","rf_cad")) {
  cols <- colnames(df)
  tmp <- data.frame(cols = df[,i])
  colnames(tmp) <- c("pval")
  df$BH <- stats::p.adjust(tmp$pval, method = "BH")
  df$evidence <- df$BH<threshold
  colnames(df) <- c(cols,paste0(i,c("_adjust","_evidence")))
}

df <- df[,c("id",
            paste0("rf_t2d",c("","_adjust","_evidence")),
            paste0("t2d_rf",c("","_adjust","_evidence")),
            paste0("rf_pad",c("","_adjust","_evidence")),
            paste0("rf_cad",c("","_adjust","_evidence")))]

df <- merge(master, df, by = "id", all.x = TRUE)

data.table::fwrite(df,"output/evidence_summary.csv")