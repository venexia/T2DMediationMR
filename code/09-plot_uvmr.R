rm(list=ls())
graphics.off()

library(magrittr)

# Specify paths ----------------------------------------------------------------

source("code/specify_paths.R", echo = TRUE)
source("code/fn-uvmr_plot.R", echo = TRUE)

# Load risk factor data --------------------------------------------------------

risk_factors <- data.table::fread("data/risk_factors_10snps.csv", data.table = FALSE)
outcomes <- data.table::fread("raw/outcomes.csv", data.table = FALSE)
gwas <- rbind(risk_factors, outcomes)

gwas$trait <- gsub("  automated.*","",gwas$trait)
gwas$trait <- ifelse(gwas$trait=="Alcohol intake frequency.","Alcohol intake frequency",gwas$trait)

binary_traits <- outcomes$id

# Load evidence summary --------------------------------------------------------

evidence <-  data.table::fread("output/evidence_summary.csv", data.table = FALSE)

# Load UVMR results ------------------------------------------------------------

results_fwd <- data.table::fread("output/results_fwd.csv", data.table = FALSE)
results_bkwd <- data.table::fread("output/results_bkwd.csv", data.table = FALSE)
df <- rbind(results_fwd, results_bkwd)
df$id.exposure <- df$exposure
df$id.outcome <- gsub(".*id:","",df$outcome)
df[,c("exposure","outcome","snps_in","snps_mr")] <- NULL
df <- df[df$id.exposure %in% gwas$id,]
df <- df[df$id.outcome %in% gwas$id,]

# Tidy -------------------------------------------------------------------------

rm(risk_factors, outcomes, results_bkwd, results_fwd)

# Annotate results -------------------------------------------------------------

trait_labels <- unique(gwas[,c("id","trait")])

colnames(trait_labels) <- c("id.exposure","exposure")
df <- merge(df, trait_labels, by = "id.exposure", all.x = TRUE)

colnames(trait_labels) <- c("id.outcome","outcome")
df <- merge(df, trait_labels, by = "id.outcome", all.x = TRUE)

# Add confidence intervals -----------------------------------------------------

df$lci <- df$b - qnorm(0.975)*df$se
df$uci <- df$b + qnorm(0.975)*df$se

# Determine whether to use beta or odds ratio ----------------------------------

df$est <- ifelse(df$id.outcome %in% binary_traits, exp(df$b), df$b)
df$est_lci <- ifelse(df$id.outcome %in% binary_traits, exp(df$lci), df$lci)
df$est_uci <- ifelse(df$id.outcome %in% binary_traits, exp(df$uci), df$uci)

# Plot for each outcome --------------------------------------------------------

df_plot <- df[df$id.outcome=="t2d" & !(df$id.exposure %in% c("cad","pad","t2d_linear")),]
df_plot <- merge(df_plot, evidence[,c("id","rf_t2d_evidence")], by.x = "id.exposure", by.y = "id", all.x = TRUE)
df_plot <- dplyr::rename(df_plot, evidence = rf_t2d_evidence)
uvmr_plot(dat = df_plot, type = "outcome", trait = "t2d")

df_plot <- df[df$id.exposure=="t2d" & !(df$id.outcome %in% c("cad","pad","t2d_linear")),]
df_plot <- merge(df_plot, evidence[,c("id","t2d_rf_evidence")], by.x = "id.outcome", by.y = "id", all.x = TRUE)
df_plot <- dplyr::rename(df_plot, evidence = t2d_rf_evidence)
uvmr_plot(dat = df_plot,type = "exposure", trait = "t2d")

df_plot <- df[df$id.outcome=="cad" & !(df$id.exposure %in% c("pad","t2d_linear","t2d")),]
df_plot <- merge(df_plot, evidence[,c("id","rf_cad_evidence")], by.x = "id.exposure", by.y = "id", all.x = TRUE)
df_plot <- dplyr::rename(df_plot, evidence = rf_cad_evidence)
uvmr_plot(dat = df_plot, type = "outcome", trait = "cad")

df_plot <- df[df$id.outcome=="pad" & !(df$id.exposure %in% c("cad","t2d_linear","t2d")),]
df_plot <- merge(df_plot, evidence[,c("id","rf_pad_evidence")], by.x = "id.exposure", by.y = "id", all.x = TRUE)
df_plot <- dplyr::rename(df_plot, evidence = rf_pad_evidence)
uvmr_plot(dat = df_plot,type = "outcome",trait = "pad")
