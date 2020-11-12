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
  df_plot$order <- seq(1:nrow(df_plot))
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

  df_plot2 <- df[df$outcome==i &
                   df$exposure %in% exposure_list,]

  df_plot2 <- merge(df_plot2,
                    df_plot[,c("exposure","outcome","order")],
                    by = c("exposure","outcome"))

  df_plot2$method <- ifelse(df_plot2$method %in% c("Wald ratio","Inverse variance weighted"),
                            "Inverse variance weighted / Wald ratio",
                            df_plot2$method)

  df_plot2 <- df_plot2[!(df_plot2$nsnp<9 & df_plot2$method %in% setdiff(unique(df_plot2$method),"Inverse variance weighted / Wald ratio")),]

  plotlabs1 <- df_plot2$exposure_long
  plotlabs2 <- df_plot2$exposure_plot
  
  plot_min <- 0.13
  plot_max <- 8
  
  df_plot2$lci_plot <- ifelse(df_plot2$lci_or<plot_min,plot_min,df_plot2$lci_or)
  df_plot2$uci_plot <- ifelse(df_plot2$uci_or>plot_max,plot_max,df_plot2$uci_or)
  df_plot2$or_plot <- ifelse(df_plot2$or<plot_min | df_plot2$or>plot_max, NA, df_plot2$or)
  
  ggplot2::ggplot(data = df_plot2,
                  mapping = ggplot2::aes(x = or_plot, y = forcats::fct_rev(method), color = method)) +
    ggplot2::geom_vline(xintercept=1, col = "dark gray") +
    ggplot2::geom_point(shape = 15, size = 0.5) +
    ggplot2::geom_linerange(ggplot2::aes(xmin = lci_plot, xmax = uci_plot), alpha = 0.6, size = 1) +
    ggplot2::scale_x_continuous(trans = "log", breaks = c(0.25,1,4),lim = c(plot_min,plot_max)) +
    ggplot2::scale_color_manual(values = c("#e41a1c","#377eb8","#4daf4a","#984ea3","#ff7f00")) +
    ggplot2::theme_minimal() +
    ggplot2::labs(x = paste0("Odds ratio and 95% confidence interval for\nthe effect of the feature on ",i_long), y = "") +
    ggplot2::theme(panel.grid.major.y = ggplot2::element_blank(),
                   panel.grid.minor.y = ggplot2::element_blank(),
                   axis.text.y = ggplot2::element_blank(),
                   axis.text.x = ggplot2::element_text(size=8),
                   text = ggplot2::element_text(size=8),
                   legend.title = ggplot2::element_blank(),
                   legend.position = "bottom") +
    ggplot2::facet_wrap(.~exposure_long, ncol = ifelse(i=="t2d",4,5))
  
  ggplot2::ggsave(filename = paste0("output/uvmr_",i,"_sa.jpeg"),
                  dpi = 300, height = 279.4, width = 215.9, unit = "mm", scale = 1)
  
}
