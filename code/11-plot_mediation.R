rm(list=ls())
graphics.off()

# Load feature info ------------------------------------------------------------

features <- data.table::fread("raw/gwas.csv",
                              stringsAsFactors = FALSE,
                              data.table = FALSE)

features <- features[features$feature==TRUE,]

# Load MVMR results ------------------------------------------------------------

df <- data.table::fread("output/all_results.csv", data.table = FALSE)

# Calculate proportion mediated ------------------------------------------------

df <- merge(df,features[,c("trait","trait_long")],by.x = "exposure",by.y = "trait")

df <- df[df$outcome %in% c("cad","pad") & df$effect %in% c("direct","total","indirect"),c("trait_long","outcome","effect","estimate","se")]

df <- tidyr::pivot_wider(df, names_from = "effect", values_from = c("estimate","se"))

df <- na.omit(df)

df <- df[abs(sign(df$estimate_total)+sign(df$estimate_indirect)+sign(df$estimate_direct))==3,]
df$proportion_mediated <- df$estimate_indirect / df$estimate_total
df$proportion_mediated_se <- sqrt((df$se_indirect)^2 + (df$se_total)^2)
df$lci <- df$proportion_mediated - qnorm(0.975)*df$proportion_mediated_se
df$uci <- df$proportion_mediated + qnorm(0.975)*df$proportion_mediated_se

# Format trait labels ----------------------------------------------------------

df$trait_long <- gsub("adjusted","adj.",df$trait_long)
df$trait_long <- gsub("distribution","dist.",df$trait_long)
df$trait_long <- gsub("concentration","con.",df$trait_long)
df$trait_long <- gsub("females","F",df$trait_long)
df$trait_long <- gsub("males","M",df$trait_long)

df$trait_plot <- paste0(df$trait_long,"\n",
                        ifelse(nchar(signif(100*df$proportion_mediated, digits = 2))<2, paste0(signif(100*df$proportion_mediated, digits = 2),".0"),signif(100*df$proportion_mediated, digits = 2)),
                        "% (95% CI: ",
                        ifelse(nchar(signif(100*df$lci, digits = 2))<2,paste0(signif(100*df$lci, digits = 2),".0"),signif(100*df$lci, digits = 2)),
                        " to ",
                        ifelse(nchar(signif(100*df$uci, digits = 2))<2,paste0(signif(100*df$uci, digits = 2),".0"),signif(100*df$uci, digits = 2)),
                        ")")

# Make plot for each outcome ---------------------------------------------------

for (i in c("pad","cad")) {
  
  i_long <- ifelse(i == "pad", 
                   "peripheral artery disease", 
                   ifelse(i == "cad", 
                          "coronary artery disease",
                          NA))
  
  ggplot2::ggplot(data = df[df$outcome==i & df$proportion_mediated<=1,], 
                  mapping = ggplot2::aes(y = forcats::fct_reorder(trait_plot, proportion_mediated), x = proportion_mediated)) +
    ggplot2::geom_linerange(ggplot2::aes(xmin = lci, xmax = uci), alpha = 0.6, position=ggplot2::position_dodge(width=0.5), colour = "dark grey") +
    ggplot2::geom_point(shape = 15, size = 0.15, position=ggplot2::position_dodge(width=0.5), colour = "dark grey") +
    ggplot2::theme_minimal() +
    ggplot2::scale_x_continuous(lim = c(0,1), breaks = seq(0,1,0.2), labels = 100*seq(0,1,0.2)) +
    ggplot2::labs(x = paste0("Proportion of the effect of\nthe feature on ",i_long,"\nmediated by type 2 diabetes"), y = "") +
    ggplot2::guides(col = ggplot2::guide_legend(ncol = 1)) +
    ggplot2::theme(panel.grid.major.y = ggplot2::element_blank(), 
                   panel.grid.minor = ggplot2::element_blank(),
                   axis.text = ggplot2::element_text(size=8),
                   text = ggplot2::element_text(size=8),
                   legend.position = "none")
  
  ggplot2::ggsave(filename = paste0("output/mvmr_",i,".tiff"),
                  dpi = 300,width = 210, height = 297, unit = "mm", scale = 0.5)
  
}