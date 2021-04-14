rm(list=ls())
graphics.off()

library(magrittr)

# Specify paths ----------------------------------------------------------------

source("code/specify_paths.R", echo = TRUE)
source("code/fn-uvmr_plot.R", echo = TRUE)

# Load feature data ------------------------------------------------------------

gwas <- data.table::fread("raw/gwas.csv",
                          stringsAsFactors = FALSE,
                          select = c("trait","trait_long"),
                          data.table = FALSE)

# Load UVMR results ------------------------------------------------------------

df <- data.table::fread("output/results.csv",
                        stringsAsFactors = FALSE,
                        data.table = FALSE)

# Annotate results -------------------------------------------------------------

df <- merge(gwas, df, by.x = "trait", by.y = "exposure", all.x = TRUE)

df <- df %>%
  dplyr::mutate(exposure = trait,
                exposure_long = trait_long) %>%
  dplyr::select(-trait, -trait_long)

df <- merge(gwas, df, by.x = "trait", by.y = "outcome", all.x = TRUE)

df <- df %>%
  dplyr::mutate(outcome = trait,
                outcome_long = trait_long) %>%
  dplyr::select(-trait, -trait_long)

# Add confidence intervals -----------------------------------------------------

df$lci <- df$b - qnorm(0.975)*df$se
df$uci <- df$b + qnorm(0.975)*df$se

# Convert to odds ratios -------------------------------------------------------

df$or <- exp(df$b)
df$lci_or <- exp(df$lci)
df$uci_or <- exp(df$uci)

# Plot for each outcome --------------------------------------------------------

uvmr_plot(dat = df[df$exposure=="t2d" & df$outcome!="t2d_linear",], 
          type = "exposure", 
          trait = "t2d")

uvmr_plot(dat = df[df$outcome=="t2d" & !(df$exposure %in% c("cad","pad","t2d_linear")),], 
          type = "outcome", 
          trait = "t2d")

uvmr_plot(dat = df[df$outcome=="cad" & !(df$exposure %in% c("pad","t2d_linear")),], 
          type = "outcome", 
          trait = "cad")

uvmr_plot(dat = df[df$outcome=="pad" & !(df$exposure %in% c("cad","t2d_linear")),],
          type = "outcome", 
          trait = "pad")

uvmr_plot(dat = df[df$exposure=="t2d_linear" & df$outcome!="t2d",], 
          type = "exposure", 
          trait = "t2d_linear")

uvmr_plot(dat = df[df$outcome=="t2d_linear" & !(df$exposure %in% c("cad","pad","t2d")),], 
          type = "outcome", 
          trait = "t2d_linear")