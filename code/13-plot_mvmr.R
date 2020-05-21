rm(list=ls())
graphics.off()

features <- data.table::fread("data/features_sources.csv",
                              stringsAsFactors = FALSE,
                              data.table = FALSE)

order <- data.table::fread("output/t2d_assoc_rank.csv",
                           stringsAsFactors = FALSE,
                           data.table = FALSE)

df <- data.table::fread("output/mvmr_results_indirect.csv",
                        stringsAsFactors = FALSE,
                        data.table = FALSE)

df <- df[df$effect %in% c("direct","total","indirect_dc"),]
df <- df[df$exposure!="t2d",]
df <- df[df$outcome!="t2d",]

df <- merge(df,features[,c("trait","trait_long")],by.x = "exposure",by.y = "trait",all.y = TRUE)
df <- merge(df,order,by = "trait_long")

df$lci <- df$estimate - qnorm(0.975)*df$se
df$uci <- df$estimate + qnorm(0.975)*df$se

df$or <- exp(df$estimate)
df$lci_or <- exp(df$lci)
df$uci_or <- exp(df$uci)

df$trait_long <- gsub("adjusted","adj.",df$trait_long)
df$trait_long <- gsub("distribution","dist.",df$trait_long)
df$trait_long <- gsub("concentration","con.",df$trait_long)
df$trait_long <- gsub("females","F",df$trait_long)
df$trait_long <- gsub("males","M",df$trait_long)

for (i in c("pad","cad")) {
  
  df_plot <- df[df$outcome==i,]
  
  ggplot2::ggplot(data = df_plot, mapping = ggplot2::aes(x = forcats::fct_reorder(trait_long, order), y = or, color = effect)) +
    ggplot2::geom_hline(yintercept=1, col = "dark gray") +
    ggplot2::geom_linerange(ggplot2::aes(ymin = lci_or, ymax = uci_or), alpha = 0.6, position=ggplot2::position_dodge(width=0.5)) +
    ggplot2::geom_linerange(ggplot2::aes(ymin = or-(1e-3), ymax = or+(1e-3)), alpha = 1, position=ggplot2::position_dodge(width=0.5)) +
    ggplot2::scale_y_continuous(trans = "log", breaks = c(0.12,0.25,0.5,1,2,4,8,16),lim = c(0.06,32)) +
    ggplot2::scale_color_manual(breaks = c("direct","indirect_dc","total"), values = c("#e41a1c","#377eb8","#984ea3"), labels = c("Effect not through type 2 diabetes","Effect through type 2 diabetes","Total effect")) +
    ggplot2::theme_minimal() +
    ggplot2::labs(y = "OR and 95% CI", x = "") +
    ggplot2::theme(panel.grid.major.x = ggplot2::element_blank(), 
                   panel.grid.minor = ggplot2::element_blank(),
                   axis.text.x = ggplot2::element_text(angle = 90,hjust = 1,vjust = 0.5),
                   axis.text = ggplot2::element_text(size=8),
                   text = ggplot2::element_text(size=8),
                   legend.position = "bottom",
                   legend.title = ggplot2::element_blank())
  
  ggplot2::ggsave(filename = paste0("output/mvmr_",i,".tiff"),
                  dpi = 300,width = 33.87, height = 19.05, unit = "cm", scale = 0.6)
  
}