rm(list=ls())
graphics.off()

# Load source data info --------------------------------------------------------

gwas <- data.table::fread("raw/gwas.csv",
                          stringsAsFactors = FALSE,
                          data.table = FALSE)

gwas <- gwas[gwas$source!="exclude_feature" & !grepl("ukb",gwas$ieugwas) & !(gwas$trait %in% c("t2d","pad","cad")),]
rownames(gwas) <- NULL

results <- data.table::fread("output/results.csv", data.table = FALSE)
results <- results[results$outcome=="t2d" & (results$exposure %in% gwas$trait),c("exposure","outcome","method","b","se","nsnp")]
results$type <- "logistic"

results_linear <- data.table::fread("output/results_linear.csv", data.table = FALSE)
results_linear <- results_linear[,c("exposure","outcome","method","b","se","nsnp")]
results_linear$type <- "linear"

results_linear_convert <- data.table::fread("output/results_linear_convert.csv", data.table = FALSE)
results_linear_convert <- results_linear_convert[,c("exposure","outcome","method","b","se","nsnp")]
results_linear_convert$type <- "linear_convert"

df <- rbind(results, results_linear_convert)

df$lci <- df$b - qnorm(0.975)*df$se
df$uci <- df$b + qnorm(0.975)*df$se

df <- df[df$method=="Inverse variance weighted",]
df$or <- exp(df$b)
df$or_lci <- exp(df$lci)
df$or_uci <- exp(df$uci)

ggplot2::ggplot(data = df, 
                mapping = ggplot2::aes(x = or, y = forcats::fct_rev(exposure), color = type)) +
  ggplot2::geom_vline(xintercept=1, col = "dark gray") +
  ggplot2::geom_linerange(ggplot2::aes(xmin = or_lci, xmax = or_uci), alpha = 0.5, position=ggplot2::position_dodge(width=0.5)) +
  ggplot2::geom_point(shape = 15, size = 0.15, position=ggplot2::position_dodge(width=0.5)) +
  ggplot2::scale_x_continuous(trans = "log", lim = c(0.13,64), breaks = c(0.13,0.25,0.5,1,2,4,8,16,32,64)) +
  ggplot2::scale_y_discrete(position = "left") +
  ggplot2::scale_color_manual(breaks = c("logistic","linear_convert"), values = c("#e41a1c","#377eb8"), labels = c("logistic","linear")) +
  ggplot2::theme_minimal() +
  ggplot2::labs(x = "Odds ratio and 95% confidence interval for the effect of the feature on type 2 diabetes", y = "", color = "Type 2 diabetes GWAS model") +
  ggplot2::guides(col = ggplot2::guide_legend(nrow = 1)) +
  ggplot2::theme(panel.grid.major.y = ggplot2::element_blank(), 
                 panel.grid.minor = ggplot2::element_blank(),
                 axis.text = ggplot2::element_text(size=8),
                 text = ggplot2::element_text(size=8),
                 strip.text = ggplot2::element_text(angle = 0),
                 legend.position = "bottom")
