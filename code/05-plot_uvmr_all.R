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

df$exposure_plot <- paste0("OR: ",
                        ifelse(nchar(signif(df$or, digits = 2))<2, paste0(signif(df$or, digits = 2),".0"),signif(df$or, digits = 2)),
                        "; 95% CI: ",
                        ifelse(nchar(signif(df$lci_or, digits = 2))<2, paste0(signif(df$lci_or, digits = 2),".0"),signif(df$lci_or, digits = 2)),
                        " to ",
                        ifelse(nchar(signif(df$uci_or, digits = 2))<2, paste0(signif(df$uci_or, digits = 2),".0"),signif(df$uci_or, digits = 2)))

# Plot type 2 diabetes results -------------------------------------------------

for (i in c("t2d","cad","pad")) {

  i_long <- ifelse(i == "pad", 
                   "peripheral artery disease", 
                   ifelse(i == "cad", 
                          "coronary artery disease",
                          ifelse(i == "t2d", 
                                 "type 2 diabates",
                                 NA)))
  
  exposure_list <- ifelse(i == "t2d", "features", "t2d_evidence")
  exposure_list <- get(exposure_list)
  
  df_plot <- df[df$method %in% c("Wald ratio","Inverse variance weighted") & 
                  df$outcome==i & 
                  df$exposure %in% exposure_list,]
  
  df_plot <- df_plot[order(df_plot$or),]
  plotlabs1 <- df_plot$exposure_long
  plotlabs2 <- df_plot$exposure_plot
  
ggplot2::ggplot(data = df_plot, 
                mapping = ggplot2::aes(x = or, y = 1:nrow(df_plot))) +
  ggplot2::geom_vline(xintercept=1, col = "dark gray") +
  ggplot2::geom_point(shape = 15, size = 0.5, col = "black") +
  ggplot2::geom_linerange(ggplot2::aes(xmin = lci_or, xmax = uci_or), alpha = 0.6, size = 1, col = "dark grey") +
  ggplot2::scale_x_continuous(trans = "log", breaks = c(0.12,0.25,0.5,1,2,4,8,16,32),lim = c(0.12,32)) +
  ggplot2::scale_y_continuous(breaks = 1:length(plotlabs1),
                              labels = plotlabs1,
                              sec.axis = ggplot2::sec_axis(~.,
                                                  breaks = 1:length(plotlabs2),
                                                  labels = plotlabs2),
                              expand = c(0,0.6)) +
  ggplot2::theme_minimal() +
  ggplot2::labs(x = paste0("Odds ratio and 95% confidence interval for\nthe effect of the feature on ",i_long), y = "") +
  ggplot2::theme(panel.grid.major.y = ggplot2::element_blank(), 
                 panel.grid.minor = ggplot2::element_blank(),
                 axis.text = ggplot2::element_text(size=8),
                 text = ggplot2::element_text(size=8))

ggplot2::ggsave(filename = paste0("output/uvmr_",i,".tiff"),
                dpi = 300,width = 210, height = 297, unit = "mm", scale = 0.8)


}
