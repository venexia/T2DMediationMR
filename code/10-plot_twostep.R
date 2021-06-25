rm(list=ls())
graphics.off()

# Load feature info ------------------------------------------------------------

gwas <- data.table::fread("data/gwas.csv", data.table = FALSE)

# Load MVMR results ------------------------------------------------------------

df <- data.table::fread("output/twostep_results.csv", data.table = FALSE)

df <- df[df$exposure!="t2d_linear",]
df <- df[df$mediator!="t2d_linear",]

df <- merge(df,gwas[,c("id","trait")],by.x = "exposure",by.y = "id", all.x = TRUE)
df <- dplyr::rename(df, exposure_long = trait)
df <- merge(df,gwas[,c("id","trait")],by.x = "mediator",by.y = "id", all.x = TRUE)
df <- dplyr::rename(df, mediator_long = trait)

df$lci <- df$estimate - qnorm(0.975)*df$se
df$uci <- df$estimate + qnorm(0.975)*df$se

df$or <- exp(df$estimate)
df$lci_or <- exp(df$lci)
df$uci_or <- exp(df$uci)

df$exposure_plot <- paste0(toupper(substr(df$effect,1,1)),substr(df$effect,2,nchar(df$effect)),", OR: ",
                           sprintf("%.2f",df$or),
                           ", 95% CI: ",
                           sprintf("%.2f",df$lci_or),
                           " to ",
                           sprintf("%.2f",df$uci_or))

df$facet_plot <- paste0("Exposure: ",df$exposure_long, "; Mediator: ", df$mediator_long)

# Make plot for each outcome ---------------------------------------------------

for (outcome in c("pad","cad")) {
  
  outcome_long <- ifelse(outcome == "pad", "peripheral artery disease", 
                   ifelse(outcome == "cad", "coronary artery disease",
                          NA))
  
  ggplot2::ggplot(data = df[df$outcome==outcome,], 
                  mapping = ggplot2::aes(y = exposure_plot, x = or, color = effect)) +
    ggplot2::geom_vline(xintercept=1, col = "dark gray") +
    ggplot2::geom_linerange(ggplot2::aes(xmin = lci_or, xmax = uci_or), alpha = 0.5, position=ggplot2::position_dodge(width=0.5)) +
    ggplot2::geom_point(shape = 15, size = 0.15, position=ggplot2::position_dodge(width=0.5)) +
    ggplot2::scale_x_continuous(trans = "log", breaks = c(0.5,1,2,4),lim = c(0.5,6)) +
    ggplot2::scale_y_discrete(position = "left") +
    ggplot2::scale_color_manual(breaks = c("direct","indirect","total"), values = c("#e41a1c","#377eb8","#984ea3"), labels = c("Independent of type 2 diabetes","Through type 2 diabetes","Total")) +
    ggplot2::theme_minimal() +
    ggplot2::labs(x = paste0("Odds ratio and 95% confidence interval for\nthe effect on ",outcome_long), y = "", color = "Effect") +
    ggplot2::guides(col = ggplot2::guide_legend(ncol = 1)) +
    ggplot2::theme(panel.grid.major.y = ggplot2::element_blank(), 
                   panel.grid.minor = ggplot2::element_blank(),
                   axis.text = ggplot2::element_text(size=8),
                   text = ggplot2::element_text(size=8),
                   strip.text = ggplot2::element_text(angle = 0),
                   legend.position = "none") +
    ggplot2::facet_wrap(facet_plot~., scales = "free_y", 
                        nrow = length(unique(df[df$outcome==outcome,]$facet_plot)), 
                        strip.position = "top")
  
  ggplot2::ggsave(filename = paste0("output/mvmr_",outcome,".jpeg"),
                  dpi = 300, width = 210, 
                  height = 297*(nrow(unique(df[df$outcome==outcome,c("exposure","mediator")]))/nrow(unique(df[df$outcome=="cad",c("exposure","mediator")]))), 
                  unit = "mm", scale = 0.7)
  
}
