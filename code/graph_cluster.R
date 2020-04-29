library(ggplot2)

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

df <- merge(clusters, features, by = "trait_udler", all.x = TRUE)

results <- data.table::fread("output/mvmr_results.csv",
                             data.table = FALSE,
                             stringsAsFactors = FALSE)

df <- merge(df, results, by.x = "trait", by.y = "feature", all.x = TRUE)

df$lci <- df$beta - qnorm(0.975)*df$se
df$uci <- df$beta + qnorm(0.975)*df$se

df$t2d <- ifelse(df$exposure=="type 2 diabetes",TRUE,FALSE)

tmp <- df[df$analysis=="MVMR" & df$t2d==FALSE,]
tmp <- tmp[!is.na(tmp$trait_long),]
  
ggplot(data = tmp, 
       mapping = aes(x = forcats::fct_rev(trait_long), y = beta)) +
  geom_point(aes(shape = direct)) + 
  geom_linerange(aes(ymin = lci, ymax = uci)) +
  geom_hline(yintercept=0, linetype = 2) +
  theme_minimal() +
  labs(y = "Beta and 95% confidence interval", x = "") +
  facet_grid(name~direct, scales = "free_y", space = "free_y") +
  coord_flip() +
  theme(panel.grid.major.y = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text=element_text(size=8),
        panel.background = element_rect(linetype = "solid"),
        text=element_text(size=8),
        legend.position="none")
