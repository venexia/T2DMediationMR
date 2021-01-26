rm(list=ls())
graphics.off()

# Specify paths ----------------------------------------------------------------

source("code/specify_paths.R")

# Load feature data ------------------------------------------------------------

gwas <- data.table::fread("raw/gwas.csv",
                          stringsAsFactors = FALSE,
                          select = c("trait","trait_long","source"),
                          data.table = FALSE)

colnames(gwas) <- c("exposure","exposure_long","exposure_source")

# Identify features for plotting -----------------------------------------------

gwas <- gwas[gwas$exposure_source %in% c("original_feature","new_feature"),]

features <- gwas$exposure

# Load UVMR results ------------------------------------------------------------

df <- data.table::fread("output/results.csv",
                        stringsAsFactors = FALSE,
                        data.table = FALSE)

df <- df[df$exclude==FALSE,]

# Annotate results -------------------------------------------------------------

df <- merge(gwas, df, by = c("exposure"),all.x = TRUE)

# Add confidence intervals -----------------------------------------------------

df$lci <- df$b - qnorm(0.975)*df$se
df$uci <- df$b + qnorm(0.975)*df$se

# Convert to odds ratios -------------------------------------------------------

df$or <- exp(df$b)
df$lci_or <- exp(df$lci)
df$uci_or <- exp(df$uci)

# Plot for each outcome --------------------------------------------------------

for (outcome in c("t2d","cad","pad")) {
  
  # Label outcome --------------------------------------------------------------
  
  outcome_long <- ""
  outcome_long <- ifelse(outcome=="t2d","type 2 diabetes",outcome_long)
  outcome_long <- ifelse(outcome=="pad","peripheral artery disease",outcome_long)
  outcome_long <- ifelse(outcome=="cad","coronary artery disease",outcome_long)
  
  # Restricy data to that needed for plotting ----------------------------------
  
  df_plot <- df[df$outcome==outcome,]
  df_plot <- df_plot[,c("exposure","exposure_long","method","or","lci_or","uci_or","pval","nsnp")]
  df_plot <- unique(df_plot)
  df_plot <- na.omit(df_plot)
  
  # Convert data to wide (i.e. one row per exposure) ---------------------------
  
  df_plot$method <- ifelse(df_plot$method %in% c("Inverse variance weighted","Wald ratio"),"main",df_plot$method)
  
  df_plot$method <- gsub(" ","_",df_plot$method)
  
  df_plot <- tidyr::pivot_wider(df_plot,
                                id_cols = c("exposure","exposure_long","nsnp"),
                                names_from = "method",
                                names_sep = "-",
                                values_from = c("or","lci_or","uci_or","pval"),
                                values_fill = NA)
  
  # Convert data to long with main result recorded alongside sensitivity -------
  
  df_plot <- tidyr::pivot_longer(df_plot,
                                 cols = c(paste0("or-",c("MR_Egger","Simple_mode","Weighted_mode","Weighted_median")),
                                          paste0("lci_or-",c("MR_Egger","Simple_mode","Weighted_mode","Weighted_median")),
                                          paste0("uci_or-",c("MR_Egger","Simple_mode","Weighted_mode","Weighted_median")),
                                          paste0("pval-",c("MR_Egger","Simple_mode","Weighted_mode","Weighted_median"))),
                                 names_to = "sensitivity")
  
  # Convert or, lci and uci data back to wide ----------------------------------
  
  df_plot$stat <- gsub("-.*","",df_plot$sensitivity)
  
  df_plot$sensitivity <- gsub(".*-","",df_plot$sensitivity)
  
  df_plot <- tidyr::pivot_wider(df_plot,
                                id_cols = c("exposure","exposure_long","nsnp","or-main","lci_or-main","uci_or-main","pval-main","sensitivity"),
                                names_from = "stat", 
                                values_from = "value")
  
  # Format data ahead of plotting ----------------------------------------------
  
  df_plot$sensitivity<- gsub("_"," ",df_plot$sensitivity)
  colnames(df_plot) <- c("exposure","exposure_long","nsnp","main_or","main_lci","main_uci","main_pval","sensitivity_analysis","sensitivity_or","sensitivity_lci","sensitivity_uci","sensitivity_pval")
  
  df_plot$exposure_plot <- paste0("OR: ",
                                  ifelse(nchar(signif(df_plot$main_or, digits = 2))<2, paste0(signif(df_plot$main_or, digits = 2),".0"),signif(df_plot$main_or, digits = 2)),
                                  "; 95% CI: ",
                                  ifelse(nchar(signif(df_plot$main_lci, digits = 2))<2, paste0(signif(df_plot$main_lci, digits = 2),".0"),signif(df_plot$main_lci, digits = 2)),
                                  " to ",
                                  ifelse(nchar(signif(df_plot$main_uci, digits = 2))<2, paste0(signif(df_plot$main_uci, digits = 2),".0"),signif(df_plot$main_uci, digits = 2)),
                                  "; ",
                                  df_plot$nsnp," SNPs")
  
  # Order main results by effect size ------------------------------------------
  
  df_plot$yorder <- rank(df_plot$main_or)
  
  # Specify breaks and labels --------------------------------------------------
  
  log_breaks <- c(0.03,0.06,0.12,0.25,0.5,1,2,4,8,16,32,64,128,256)
  log_labels <- c("0.03","0.06","0.12","0.25","0.5","1","2","4","8","16","32","64","128","256")
  
  # Plot main results ----------------------------------------------------------
  
  ggplot2::ggplot(data = df_plot,
                  mapping = ggplot2::aes(x = main_or, y = yorder)) +
    ggplot2::geom_vline(xintercept=1, col = "dark grey") +
    ggplot2::geom_linerange(ggplot2::aes(xmin = main_lci, xmax = main_uci), alpha = 0.1, size = 1, col = "#084594") +
    ggplot2::geom_point(shape = 15, size = 0.5, col = "#084594") +
    ggplot2::scale_x_continuous(trans = "log", 
                                breaks = log_breaks,
                                labels = log_labels,
                                lim = c(round(min(df_plot$main_lci),2)-0.01,round(max(df_plot$main_uci),2)+1)) +
    ggplot2::scale_y_continuous(breaks = df_plot$yorder,
                                labels = df_plot$exposure_long,
                                sec.axis = ggplot2::sec_axis(~.,
                                                             breaks = df_plot$yorder,
                                                             labels = df_plot$exposure_plot),
                                expand = c(0,0.6)) +
    ggplot2::theme_minimal() +
    ggplot2::labs(x = paste0("Odds ratio and 95% confidence interval for\nthe effect of the feature on ",outcome_long), y = "") +
    ggplot2::theme(panel.grid.major.y = ggplot2::element_blank(),
                   panel.grid.minor = ggplot2::element_blank(),
                   axis.text = ggplot2::element_text(size=8),
                   text = ggplot2::element_text(size=8))
  
  ggplot2::ggsave(filename = paste0("output/uvmr_",outcome,".jpeg"),
                  dpi = 300, width = 210, height = 297, unit = "mm", scale = 0.9)
  
  # Plot sensitivity analysis results ------------------------------------------

  ggplot2::ggplot(data = df_plot[df_plot$nsnp>9,],
                  mapping = ggplot2::aes(x = main_or, y = sensitivity_or)) +
    ggplot2::geom_hline(yintercept = 1, col = "dark gray") +
    ggplot2::geom_vline(xintercept = 1, col = "dark gray") +
    ggplot2::geom_abline(intercept = 0, slope = 1, col = "dark gray") +
    ggplot2::geom_errorbarh(ggplot2::aes(xmin = main_lci, xmax = main_uci), alpha = 0.5, col = "#084594") +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = sensitivity_lci, ymax = sensitivity_uci), alpha = 0.5, col = "#084594") +
    ggplot2::scale_x_continuous(trans = "log", 
                                breaks = log_breaks,
                                labels = log_labels,
                                lim = c(round(min(df_plot[df_plot$nsnp>9,]$main_lci,df_plot[df_plot$nsnp>9,]$sensitivity_lci),2)-0.01,round(max(df_plot[df_plot$nsnp>9,]$main_uci,df_plot[df_plot$nsnp>9,]$sensitivity_uci),2)+1)) +
    ggplot2::scale_y_continuous(trans = "log", 
                                breaks = log_breaks,
                                labels = log_labels,
                                lim = c(round(min(df_plot[df_plot$nsnp>9,]$main_lci,df_plot[df_plot$nsnp>9,]$sensitivity_lci),2)-0.01,round(max(df_plot[df_plot$nsnp>9,]$main_uci,df_plot[df_plot$nsnp>9,]$sensitivity_uci),2)+1)) +
    ggplot2::theme_minimal() +
    ggplot2::labs(x = paste0("Odds ratio and 95% confidence interval for the effect of the feature on\n",outcome_long," using the inverse variance weighted method"),
                  y = paste0("Odds ratio and 95% confidence interval for the effect of the feature on\n",outcome_long," using the specified sensitivity analysis method")) +
    ggplot2::theme(panel.grid.major = ggplot2::element_blank(),
                   panel.grid.minor = ggplot2::element_blank(),
                   axis.text = ggplot2::element_text(size=8),
                   text = ggplot2::element_text(size=8),
                   legend.title = ggplot2::element_blank(),
                   legend.position = "bottom",
                   strip.placement = "outcome") +
    ggplot2::facet_wrap(.~sensitivity_analysis, strip.position = "left", scales = "free")

  ggplot2::ggsave(filename = paste0("output/uvmr_",outcome,"_sa.jpeg"),
                  dpi = 300, height = 210, width = 210, unit = "mm", scale = 1)
  
}