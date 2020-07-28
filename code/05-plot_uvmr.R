rm(list=ls())
graphics.off()

# Specify paths ----------------------------------------------------------------

source("code/specify_paths.R")

# Load feature data ------------------------------------------------------------

gwas <- data.table::fread("raw/gwas.csv",
                          stringsAsFactors = FALSE,
                          select = c("trait","trait_long","feature"),
                          data.table = FALSE)

colnames(gwas) <- c("exposure","exposure_long","feature")

features <- gwas[gwas$feature==TRUE,]$exposure

# Load UVMR results -------------------------------------------------------------

df <- data.table::fread("output/uvmr_results.csv",
                        stringsAsFactors = FALSE,
                        data.table = FALSE)

# Remove results where the exposure was not a feature --------------------------

df <- df[df$exposure %in% features,]

# Annotate results -------------------------------------------------------------

df <- merge(gwas[gwas$feature==TRUE,], df, by = c("exposure"),all.x = TRUE)

# Add confidence intervals -----------------------------------------------------

df$lci <- df$b - qnorm(0.975)*df$se
df$uci <- df$b + qnorm(0.975)*df$se

# Mark estimates in MVMR -------------------------------------------------------

df$evidence <- FALSE
df$evidence <- df$pval < 0.05
df$evidence_plot <- ifelse(df$evidence==TRUE,df$outcome,"none")

t2d_evidence <- df[df$outcome=="t2d" & df$method %in% c("Wald ratio","Inverse variance weighted") & df$evidence==TRUE,]$exposure

tmp <- df[df$method %in% c("Wald ratio","Inverse variance weighted"), c("exposure","outcome","evidence")]
tmp <- tidyr::pivot_wider(tmp, names_from = "outcome", values_from = "evidence")
tmp$mvmr <- "none"
tmp$mvmr <- ifelse(tmp$t2d==TRUE & tmp$pad==TRUE & tmp$cad==TRUE, "both", tmp$mvmr)
tmp$mvmr <- ifelse(tmp$t2d==TRUE & tmp$pad==TRUE & tmp$cad==FALSE, "pad", tmp$mvmr)
tmp$mvmr <- ifelse(tmp$t2d==TRUE & tmp$pad==FALSE & tmp$cad==TRUE, "cad", tmp$mvmr)
tmp$mvmr <- factor(tmp$mvmr)
levels(tmp$mvmr) <- c(levels(tmp$mvmr), setdiff(c("both","pad","cad","none"),levels(tmp$mvmr)))

df <- merge(df, tmp[,c("exposure","mvmr")], by = "exposure", all = TRUE)

# Rank estimates by size -------------------------------------------------------

tmp <- df[df$outcome=="t2d" & df$method %in% c("Wald ratio","Inverse variance weighted") & df$exposure %in% features,c("exposure","b")]
tmp$order <- rank(tmp$b,na.last = TRUE)
tmp$b <- NULL
data.table::fwrite(tmp,"output/t2d_assoc_rank.csv",row.names = FALSE)
df <- merge(df, tmp, by = c("exposure"), all.x = TRUE)

# Convert to odds ratios -------------------------------------------------------

df$or <- exp(df$b)
df$lci_or <- exp(df$lci)
df$uci_or <- exp(df$uci)

# Format exposure labels -------------------------------------------------------

df$exposure_long <- gsub("adjusted","adj.",df$exposure_long)
df$exposure_long <- gsub("distribution","dist.",df$exposure_long)
df$exposure_long <- gsub("concentration","con.",df$exposure_long)
df$exposure_long <- gsub("females","F",df$exposure_long)
df$exposure_long <- gsub("males","M",df$exposure_long)
df$exposure_long <- factor(df$exposure_long)

# Format outcome labels --------------------------------------------------------

df$outcome <- ifelse(df$outcome=="t2d","Type 2 diabetes",df$outcome)
df$outcome <- ifelse(df$outcome=="cad","Coronary artery disease",df$outcome)
df$outcome <- ifelse(df$outcome=="pad","Peripheral artery disease",df$outcome)
df$outcome <- factor(df$outcome)
df$outcome <- factor(df$outcome,levels(df$outcome)[c(3,1,2)])

# Plot type 2 diabetes results -------------------------------------------------

df$mvmr_t2d <- ifelse(df$exposure %in% t2d_evidence,"mvmr","none")

ggplot2::ggplot(data = df[df$method %in% c("Wald ratio","Inverse variance weighted") & df$outcome=="Type 2 diabetes",], 
                mapping = ggplot2::aes(y = forcats::fct_reorder(exposure_long, order, .desc = FALSE), 
                                       x = or)) +
  ggplot2::geom_vline(xintercept=1, col = "dark gray") +
  ggplot2::geom_point(shape = 15, size = 0.5, col = "black") +
  ggplot2::geom_linerange(ggplot2::aes(xmin = lci_or, xmax = uci_or), alpha = 0.6, size = 1, col = "dark grey") +
  ggplot2::scale_x_continuous(trans = "log", breaks = c(0.12,0.25,0.5,1,2,4,8,16,32),lim = c(0.12,32)) +
  ggplot2::theme_minimal() +
  ggplot2::labs(x = "Odds ratio and 95% confidence interval for\nthe effect of the feature on type 2 diabetes", y = "") +
  ggplot2::theme(panel.grid.major.y = ggplot2::element_blank(), 
                 panel.grid.minor = ggplot2::element_blank(),
                 axis.text = ggplot2::element_text(size=8),
                 text = ggplot2::element_text(size=8))

ggplot2::ggsave(filename = "output/uvmr_t2d.tiff",
                dpi = 300,width = 210, height = 297, unit = "mm", scale = 0.7)

# Plot UVMR results ------------------------------------------------------------

ggplot2::ggplot(data = df[df$method %in% c("Wald ratio","Inverse variance weighted") & df$exposure %in% t2d_evidence,], 
                mapping = ggplot2::aes(y = forcats::fct_reorder(exposure_long, order, .desc = FALSE), 
                                       x = or, 
                                       colour = mvmr)) +
  ggplot2::geom_vline(xintercept=1, col = "dark gray") +
  ggplot2::geom_point(shape = 15, size = 0.5) +
  ggplot2::geom_linerange(ggplot2::aes(xmin = lci_or, xmax = uci_or), alpha = 0.6, size = 1) +
  ggplot2::scale_x_continuous(trans = "log", breaks = c(0.12,0.25,0.5,1,2,4,8,16), limits = c(0.18,10)) +
  ggplot2::scale_color_manual(values = c("#984ea3","#e41a1c","#377eb8","dark grey"),
                              breaks = c("both","cad","pad","none"),
                              labels = c("MVMR for CAD and MVMR for PAD",
                                         "MVMR for CAD only",
                                         "MVMR for PAD only",
                                         "None"),
                              drop = FALSE) +
  ggplot2::theme_minimal() +
  ggplot2::labs(x = "Odds ratio and 95% confidence interval for\nthe effect of the feature on each outcome", y = "", colour = "Further analyses") +
  ggplot2::theme(panel.grid.major.y = ggplot2::element_blank(), 
                 panel.grid.minor = ggplot2::element_blank(),
                 axis.text = ggplot2::element_text(size=8),
                 text = ggplot2::element_text(size=8),
                 legend.title = ggplot2::element_text(size=8),
                 legend.position = "bottom") +
  ggplot2::facet_wrap(.~outcome)

ggplot2::ggsave(filename = "output/uvmr_results.tiff",
                dpi = 300,width = 297, height = 210, unit = "mm", scale = 0.7)