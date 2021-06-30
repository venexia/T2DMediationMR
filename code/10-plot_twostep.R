rm(list=ls())
graphics.off()

# Load feature info ------------------------------------------------------------

risk_factors <- data.table::fread("data/risk_factors_10snps.csv", data.table = FALSE)
outcomes <- data.table::fread("raw/outcomes.csv", data.table = FALSE)
gwas <- rbind(risk_factors, outcomes)
gwas$trait <- gsub("  automated.*","",gwas$trait)
gwas$trait <- ifelse(gwas$trait=="Alcohol intake frequency.","Alcohol intake frequency",gwas$trait)

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

#df$facet_plot <- paste0("Exposure: ",df$exposure_long, "; Mediator: ", df$mediator_long)

df$t2d_mediator <- df$mediator=="t2d"
df$trait <- ifelse(df$mediator=="t2d",df$exposure_long,df$mediator_long)

# Make plot for each outcome ---------------------------------------------------

for (outcome in c("pad","cad")) {
  
  outcome_long <- ifelse(outcome == "pad", "peripheral artery disease", 
                   ifelse(outcome == "cad", "coronary artery disease",
                          NA))
  
  ggplot2::ggplot(data = df[df$outcome==outcome & df$t2d_mediator==TRUE,], 
                  mapping = ggplot2::aes(y = exposure_plot, x = or, color = effect)) +
    ggplot2::geom_vline(xintercept=1, col = "dark gray") +
    ggplot2::geom_linerange(ggplot2::aes(xmin = lci_or, xmax = uci_or), alpha = 0.5, position=ggplot2::position_dodge(width=0.5)) +
    ggplot2::geom_point(shape = 15, size = 0.15, position=ggplot2::position_dodge(width=0.5)) +
    ggplot2::scale_x_continuous(trans = "log", breaks = c(0.5,1,2,4),lim = c(0.5,6)) +
    ggplot2::scale_y_discrete(position = "left") +
    ggplot2::scale_color_manual(breaks = c("direct","indirect","total"), values = c("#e41a1c","#377eb8","#984ea3"), labels = c("Independent of type 2 diabetes","Through type 2 diabetes","Total")) +
    ggplot2::theme_minimal() +
    ggplot2::labs(x = paste0("Odds ratio and 95% confidence interval for the effect of\nthe risk factor on liability to ",outcome_long,",\nmediated by liability to type 2 diabetes"), y = "", color = "Effect") +
    ggplot2::guides(col = ggplot2::guide_legend(ncol = 1)) +
    ggplot2::theme(panel.grid.major.y = ggplot2::element_blank(), 
                   panel.grid.minor = ggplot2::element_blank(),
                   axis.text = ggplot2::element_text(size=8),
                   text = ggplot2::element_text(size=8),
                   strip.text = ggplot2::element_text(angle = 0),
                   legend.position = "none") +
    ggplot2::facet_wrap(trait~., scales = "free_y", 
                        nrow = length(unique(df[df$outcome==outcome & df$t2d_mediator==TRUE,]$trait)), 
                        strip.position = "top")
  
  ggplot2::ggsave(filename = paste0("output/mvmr_",outcome,"_t2dmediator.jpeg"),
                  dpi = 300, width = 210, 
                  height = 297*((length(unique(df[df$outcome==outcome & df$t2d_mediator==TRUE,]$trait)))/8), 
                  unit = "mm", scale = 0.7)
  
  ggplot2::ggplot(data = df[df$outcome==outcome & df$t2d_mediator==FALSE,], 
                  mapping = ggplot2::aes(y = exposure_plot, x = or, color = effect)) +
    ggplot2::geom_vline(xintercept=1, col = "dark gray") +
    ggplot2::geom_linerange(ggplot2::aes(xmin = lci_or, xmax = uci_or), alpha = 0.5, position=ggplot2::position_dodge(width=0.5)) +
    ggplot2::geom_point(shape = 15, size = 0.15, position=ggplot2::position_dodge(width=0.5)) +
    ggplot2::scale_x_continuous(trans = "log", breaks = c(0.5,1,2,4),lim = c(0.5,6)) +
    ggplot2::scale_y_discrete(position = "left") +
    ggplot2::scale_color_manual(breaks = c("direct","indirect","total"), values = c("#e41a1c","#377eb8","#984ea3"), labels = c("Independent of type 2 diabetes","Through type 2 diabetes","Total")) +
    ggplot2::theme_minimal() +
    ggplot2::labs(x = paste0("Odds ratio and 95% confidence interval for the effect of\nliability to type 2 diabetes on liability to ",outcome_long,",\nmediated by the risk factor"), y = "", color = "Effect") +
    ggplot2::guides(col = ggplot2::guide_legend(ncol = 1)) +
    ggplot2::theme(panel.grid.major.y = ggplot2::element_blank(), 
                   panel.grid.minor = ggplot2::element_blank(),
                   axis.text = ggplot2::element_text(size=8),
                   text = ggplot2::element_text(size=8),
                   strip.text = ggplot2::element_text(angle = 0),
                   legend.position = "none") +
    ggplot2::facet_wrap(trait~., scales = "free_y", 
                        nrow = length(unique(df[df$outcome==outcome & df$t2d_mediator==FALSE,]$trait)), 
                        strip.position = "top")
  
  ggplot2::ggsave(filename = paste0("output/mvmr_",outcome,"_t2dexposure.jpeg"),
                  dpi = 300, width = 210, 
                  height = 297*((length(unique(df[df$outcome==outcome & df$t2d_mediator==FALSE,]$trait)))/8), 
                  unit = "mm", scale = 0.7)
  
}
