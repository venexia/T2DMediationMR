rm(list=ls())
graphics.off()

library(magrittr)

# Specify paths ----------------------------------------------------------------

source("code/specify_paths.R", echo = TRUE)
source("code/fn-uvmr_plot.R", echo = TRUE)

# Load risk factor data --------------------------------------------------------

risk_factors <- data.table::fread("data/risk_factors.csv", data.table = FALSE)
risk_factors <- data.table::fread("data/risk_factors.csv", data.table = FALSE)
outcomes <- data.table::fread("raw/outcomes.csv", data.table = FALSE)
gwas <- rbind(risk_factors, outcomes)

binary_traits <- outcomes$id

# Load UVMR results ------------------------------------------------------------

results_fwd <- data.table::fread("output/results_fwd.csv", data.table = FALSE)
results_bkwd <- data.table::fread("output/results_bkwd.csv", data.table = FALSE)
df <- rbind(results_fwd, results_bkwd)
df$id.exposure <- df$exposure
df$id.outcome <- gsub(".*id:","",df$outcome)
df[,c("exposure","outcome","snps_in","snps_mr")] <- NULL
df <- df[df$id.exposure!=df$id.outcome,]

df <- df[!(df$id.exposure=="ieu-a-1041" | df$id.outcome=="ieu-a-1041"),]

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

uvmr_plot(dat = df[df$id.outcome=="t2d" & !(df$id.exposure %in% c("cad","pad","t2d_linear")),], 
          type = "outcome", 
          trait = "t2d")

uvmr_plot(dat = df[df$id.exposure=="t2d" & !(df$id.outcome %in% c("cad","pad","t2d_linear","ieu-a-1017")),], 
          type = "exposure", 
          trait = "t2d")

uvmr_plot(dat = df[df$id.outcome=="cad" & !(df$id.exposure %in% c("pad","t2d_linear","t2d")),], 
          type = "outcome", 
          trait = "cad")

uvmr_plot(dat = df[df$id.outcome=="pad" & !(df$id.exposure %in% c("cad","t2d_linear","t2d")),],
          type = "outcome", 
          trait = "pad")
