rm(list=ls())
graphics.off()

clusters <- readxl::read_xlsx("raw/L2EU.H.mat.9_T2D_71traits_259snps_samplesizefilter_v4.xlsx",sheet = "Proinsulin",range = "M1:V143")
clusters <- tidyr::pivot_longer(clusters,cols = 2:10,)
clusters <- clusters[clusters$value>0.75,]
clusters$trait_udler <- clusters$feature
clusters$trait_udler <- gsub(".ZN_pos","",clusters$trait_udler)
clusters$trait_udler <- gsub(".ZN_neg","",clusters$trait_udler)
clusters$trait_udler <- gsub("._irnt_both","",clusters$trait_udler)
clusters$feature <- NULL

features <- readxl::read_xlsx("raw/feature_sources_ieugwas.xlsx",sheet = "All")
features <- features[,c("trait","trait_long","trait_udler")]
features$trait_udler <- gsub("UKBB.","",features$trait_udler)

clusters <- merge(clusters, features, by = "trait_udler", all.x = TRUE)

features <- data.table::fread("data/features_sources.csv",
                              stringsAsFactors = FALSE,
                              data.table = FALSE)

order <- data.table::fread("output/t2d_assoc_rank.csv",
                           stringsAsFactors = FALSE,
                           data.table = FALSE)

df_start <- data.table::fread("output/mvmr_results_indirect.csv",
                        stringsAsFactors = FALSE,
                        data.table = FALSE)

df_start <- df_start[df_start$effect %in% c("direct","total","indirect_dc"),]
df_start <- df_start[df_start$exposure!="t2d",]
df_start <- df_start[df_start$outcome!="t2d",]

for (j in c("pad","cad")) {
  
  df_outcome <- df_start[df_start$outcome==j,]
  
  for (i in 1:9) {
    
    df <- df_outcome[df_outcome$exposure %in% clusters[clusters$name==paste0("W",i),]$trait,]
    
    df <- merge(df,features[,c("trait","trait_long")],by.x = "exposure",by.y = "trait")
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
    
    ggplot2::ggplot(data = df, mapping = ggplot2::aes(x = forcats::fct_reorder(trait_long, order), y = or, color = effect)) +
      ggplot2::geom_hline(yintercept=1, col = "dark gray") +
      ggplot2::geom_linerange(ggplot2::aes(ymin = lci_or, ymax = uci_or), alpha = 0.6, size = 1, position=ggplot2::position_dodge(width=0.5)) +
      ggplot2::geom_linerange(ggplot2::aes(ymin = or-(1e-3), ymax = or+(1e-3)), alpha = 1, size = 1, position=ggplot2::position_dodge(width=0.5)) +
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
    
    ggplot2::ggsave(filename = paste0("output/mvmr_W",i,"_",j,".tiff"),
                    dpi = 300,width = 33.87, height = 19.05, unit = "cm", scale = 0.6)
    
  }
}
