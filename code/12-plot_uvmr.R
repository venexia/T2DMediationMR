rm(list=ls())
graphics.off()

features <- data.table::fread("data/features_sources.csv",
                              stringsAsFactors = FALSE,
                              data.table = FALSE)

df <- data.table::fread("output/uvmr_t2d.csv",
                                stringsAsFactors = FALSE,
                                data.table = FALSE)

df <- merge(df,features[,c("trait","trait_long")],by.x = "exposure",by.y = "trait",all.y = TRUE)

df$lci <- df$estimate - qnorm(0.975)*df$se
df$uci <- df$estimate + qnorm(0.975)*df$se

df$or <- exp(df$estimate)
df$lci_or <- exp(df$lci)
df$uci_or <- exp(df$uci)

df$or_pos <- ifelse(sign(df$estimate)==1,df$or,1/abs(df$or))
df$lci_or_pos <- ifelse(sign(df$estimate)==1,df$lci_or,1/abs(df$lci_or))
df$uci_or_pos <- ifelse(sign(df$estimate)==1,df$uci_or,1/abs(df$uci_or))

df$direction <- ""
df$direction <- ifelse(sign(df$estimate)==1,"Increases risk of type 2 diabetes",df$direction)
df$direction <- ifelse(sign(df$estimate)==-1,"Decreases risk of type 2 diabetes",df$direction)

df$excl_null <- FALSE
df$excl_null <- ifelse(df$lci_or_pos>1 & df$uci_or_pos>1,TRUE,df$excl_null)
df$excl_null <- ifelse(df$lci_or_pos<1 & df$uci_or_pos<1,TRUE,df$excl_null)

df$order <- rank(-df$or_pos,na.last = TRUE)
data.table::fwrite(df[,c("trait_long","order")],"output/t2d_assoc_rank.csv",row.names = FALSE)

df$trait_long <- gsub("adjusted","adj.",df$trait_long)
df$trait_long <- gsub("distribution","dist.",df$trait_long)
df$trait_long <- gsub("concentration","con.",df$trait_long)
df$trait_long <- gsub("females","F",df$trait_long)
df$trait_long <- gsub("males","M",df$trait_long)
df$trait_long <- factor(df$trait_long)

ggplot2::ggplot(data = df, mapping = ggplot2::aes(x = forcats::fct_reorder(trait_long, or_pos, .desc = TRUE), y = or_pos, color = direction)) +
  ggplot2::geom_hline(yintercept=1, col = "dark gray") +
  ggplot2::geom_linerange(ggplot2::aes(ymin = lci_or_pos, ymax = uci_or_pos), alpha = 0.6, size = 1, position=ggplot2::position_dodge(width=0.5)) +
  ggplot2::geom_linerange(ggplot2::aes(ymin = or_pos-(1e-3), ymax = or_pos+(1e-3)), alpha = 1, size = 1, position=ggplot2::position_dodge(width=0.5)) +
  ggplot2::scale_y_continuous(trans = "log", breaks = c(0.12,0.25,0.5,1,2,4,8,16),lim = c(0.06,32)) +
  ggplot2::scale_color_manual(breaks = c("Decreases risk of type 2 diabetes","Increases risk of type 2 diabetes"), values=c("#377eb8","#e41a1c")) +
  ggplot2::theme_minimal() +
  ggplot2::labs(y = "OR and 95% CI", x = "") +
  ggplot2::theme(panel.grid.major.x = ggplot2::element_blank(), 
                 panel.grid.minor = ggplot2::element_blank(),
                 axis.text.x = ggplot2::element_text(angle = 90,hjust = 1,vjust = 0.5),
                 axis.text = ggplot2::element_text(size=8),
                 text = ggplot2::element_text(size=8),
                 legend.position = "bottom",
                 legend.title = ggplot2::element_blank())

ggplot2::ggsave(filename = "output/uvmr_t2d.tiff",
       dpi = 300,width = 33.87, height = 19.05, unit = "cm", scale = 0.6)