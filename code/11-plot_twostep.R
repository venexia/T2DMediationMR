rm(list=ls())
graphics.off()

# Load feature info ------------------------------------------------------------

features <- data.table::fread("raw/gwas.csv",
                              stringsAsFactors = FALSE,
                              data.table = FALSE)

# Load MVMR results ------------------------------------------------------------

df <- data.table::fread("output/twostep_results.csv", data.table = FALSE)

df <- merge(df,features[,c("trait","trait_long")],by.x = "exposure",by.y = "trait")

df$lci <- df$estimate - qnorm(0.975)*df$se
df$uci <- df$estimate + qnorm(0.975)*df$se

df$or <- exp(df$estimate)
df$lci_or <- exp(df$lci)
df$uci_or <- exp(df$uci)

df$exposure_plot <- paste0(toupper(substr(df$effect,1,1)),substr(df$effect,2,nchar(df$effect)),", OR: ",
                            ifelse(nchar(signif(df$or, digits = 2))<2, paste0(signif(df$or, digits = 2),".0"),signif(df$or, digits = 2)),
                            ", 95% CI: ",
                            ifelse(nchar(signif(df$lci_or, digits = 2))<2, paste0(signif(df$lci_or, digits = 2),".0"),signif(df$lci_or, digits = 2)),
                            " to ",
                            ifelse(nchar(signif(df$uci_or, digits = 2))<2, paste0(signif(df$uci_or, digits = 2),".0"),signif(df$uci_or, digits = 2)))

# Make plot for each outcome ---------------------------------------------------

for (i in c("pad","cad")) {
  
  i_long <- ifelse(i == "pad", 
                   "peripheral artery disease", 
                   ifelse(i == "cad", 
                          "coronary artery disease",
                          NA))
  
  mvmr_features <- unique(df[df$effect=="direct" & df$outcome==i,]$exposure)
  
  df_plot <- df[df$exposure %in% mvmr_features & 
                  df$effect %in% c("direct","total","indirect") &
                  df$outcome==i,]
  
  ggplot2::ggplot(data = df_plot, 
                  mapping = ggplot2::aes(y = exposure_plot, x = or, color = effect)) +
    ggplot2::geom_vline(xintercept=1, col = "dark gray") +
    ggplot2::geom_linerange(ggplot2::aes(xmin = lci_or, xmax = uci_or), alpha = 0.5, position=ggplot2::position_dodge(width=0.5)) +
    ggplot2::geom_point(shape = 15, size = 0.15, position=ggplot2::position_dodge(width=0.5)) +
    # ggplot2::scale_x_continuous(trans = "log", breaks = c(0.12,0.25,0.5,1,2,4,8,16),lim = c(0.25,4)) +
    ggplot2::scale_x_continuous(trans = "log", breaks = c(0.35,0.75,0.5,1,1.33,2,2.86),lim = c(0.35,2.86)) +
    ggplot2::scale_y_discrete(position = "left") +
    ggplot2::scale_color_manual(breaks = c("direct","indirect","total"), values = c("#e41a1c","#377eb8","#984ea3"), labels = c("Independent of type 2 diabetes","Through type 2 diabetes","Total")) +
    ggplot2::theme_minimal() +
    ggplot2::labs(x = paste0("Odds ratio and 95% confidence interval for\nthe effect of the feature on ",i_long), y = "", color = "Effect") +
    ggplot2::guides(col = ggplot2::guide_legend(ncol = 1)) +
    ggplot2::theme(panel.grid.major.y = ggplot2::element_blank(), 
                   panel.grid.minor = ggplot2::element_blank(),
                   axis.text = ggplot2::element_text(size=8),
                   text = ggplot2::element_text(size=8),
                   strip.text = ggplot2::element_text(angle = 0),
                   legend.position = "none") +
    ggplot2::facet_wrap(trait_long~., scales = "free_y", 
                        nrow = length(unique(df_plot$exposure)), 
                        strip.position = "top")
  
  ggplot2::ggsave(filename = paste0("output/mvmr_",i,".jpeg"),
                  dpi = 300,width = 210, height = 148.5, unit = "mm", scale = 0.7)
  
}
