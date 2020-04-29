library(ggplot2)

features <- readxl::read_xlsx("raw/feature_sources_ieugwas.xlsx",sheet = "All")
features <- features[,c("trait","trait_long")]

results <- data.table::fread("output/mvmr_results.csv",
                        data.table = FALSE,
                        stringsAsFactors = FALSE)

df <- merge(features, results, by.x = "trait", by.y = "feature")

df$lci <- df$beta - qnorm(0.975)*df$se
df$uci <- df$beta + qnorm(0.975)*df$se

df$t2d <- ifelse(df$exposure=="type 2 diabetes",TRUE,FALSE)


ggplot(data = df[df$analysis=="MVMR" & df$t2d==FALSE,], 
       mapping = aes(x = forcats::fct_rev(trait_long), y = beta)) +
  geom_point() + 
  geom_linerange(aes(ymin = lci, ymax = uci)) +
  geom_hline(yintercept=0, linetype = 2) +
  theme_minimal() +
  labs(y = "Beta and 95% confidence interval", x = "") +
  facet_wrap(.~direct) +
  coord_flip() +
  theme(panel.grid.major.y = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text=element_text(size=8),
        text=element_text(size=8),
        strip.text.y = element_text(size=8,hjust = 1, angle = 180),
        legend.position="none")


ggplot(data = df[df$analysis=="MVMR" & df$t2d==TRUE,], 
       mapping = aes(x = forcats::fct_rev(trait_long), y = beta)) +
  geom_point() + 
  geom_linerange(aes(ymin = lci, ymax = uci)) +
  geom_hline(yintercept=0, linetype = 2) +
  theme_minimal() +
  labs(y = "Beta and 95% confidence interval", x = "") +
  facet_wrap(.~direct) +
  coord_flip() +
  theme(panel.grid.major.y = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text=element_text(size=8),
        text=element_text(size=8),
        strip.text.y = element_text(size=8,hjust = 1, angle = 180),
        legend.position="none")
