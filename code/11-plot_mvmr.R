rm(list=ls())
graphics.off()

# Load feature info ------------------------------------------------------------

features <- data.table::fread("raw/gwas.csv",
                              stringsAsFactors = FALSE,
                              data.table = FALSE)

features <- features[features$feature==TRUE,]

# Load UVMR results to determine which MVMR results are to be kept -------------

uvmr <- data.table::fread("output/uvmr_results.csv", data.table = FALSE)

uvmr$lci <- uvmr$b - qnorm(0.975)*uvmr$se
uvmr$uci <- uvmr$b + qnorm(0.975)*uvmr$se

uvmr <- uvmr[uvmr$method %in% c("Inverse variance weighted","Wald ratio") &
               !(uvmr$exposure %in% c("t2d","t2d_udler","cad","pad")),
             c("exposure","outcome","b","se","lci","uci")]

colnames(uvmr) <- c("exposure","outcome","estimate","se","lci","uci")

uvmr$excl_null <- sign(uvmr$lci)==sign(uvmr$uci)

uvmr[,c("lci","uci")] <- NULL

uvmr <- tidyr::pivot_wider(uvmr, 
                           names_from = "outcome", 
                           values_from = c("estimate","se","excl_null"))

cad_features <- uvmr[uvmr$excl_null_t2d==TRUE & uvmr$excl_null_cad==TRUE,]$exposure
pad_features <- uvmr[uvmr$excl_null_t2d==TRUE & uvmr$excl_null_pad==TRUE,]$exposure

# Load MVMR results ------------------------------------------------------------

mvmr <- data.table::fread("output/mvmr_results_indirect.csv", data.table = FALSE)

mvmr <- merge(mvmr,features[,c("trait","trait_long")],by.x = "exposure",by.y = "trait")

mvmr$lci <- mvmr$estimate - qnorm(0.975)*mvmr$se
mvmr$uci <- mvmr$estimate + qnorm(0.975)*mvmr$se

mvmr$or <- exp(mvmr$estimate)
mvmr$lci_or <- exp(mvmr$lci)
mvmr$uci_or <- exp(mvmr$uci)

mvmr$trait_long <- gsub("adjusted","adj.",mvmr$trait_long)
mvmr$trait_long <- gsub("distribution","dist.",mvmr$trait_long)
mvmr$trait_long <- gsub("concentration","con.",mvmr$trait_long)
mvmr$trait_long <- gsub("females","F",mvmr$trait_long)
mvmr$trait_long <- gsub("males","M",mvmr$trait_long)

# Make plot for each outcome ---------------------------------------------------

for (i in c("pad","cad")) {
  
  mvmr_features <- get(paste0(i,"_features"))
  
  df <- mvmr[mvmr$exposure %in% mvmr_features & 
               mvmr$effect %in% c("direct","total_restricted","indirect_MVMR") &
               mvmr$outcome==i,]
  
  ggplot2::ggplot(data = df, mapping = ggplot2::aes(y = forcats::fct_rev(trait_long), x = or, color = effect)) +
    ggplot2::geom_vline(xintercept=1, col = "dark gray") +
    ggplot2::geom_linerange(ggplot2::aes(xmin = lci_or, xmax = uci_or), alpha = 0.6, position=ggplot2::position_dodge(width=0.5)) +
    ggplot2::geom_linerange(ggplot2::aes(xmin = or-(1e-3), xmax = or+(1e-3)), alpha = 1, position=ggplot2::position_dodge(width=0.5)) +
    ggplot2::scale_x_continuous(trans = "log", breaks = c(0.12,0.25,0.5,1,2,4,8,16),lim = c(0.25,4)) +
    ggplot2::scale_color_manual(breaks = c("direct","indirect_MVMR","total_restricted"), values = c("#e41a1c","#377eb8","#984ea3"), labels = c("Effect not through type 2 diabetes","Effect through type 2 diabetes","Total effect")) +
    ggplot2::theme_minimal() +
    ggplot2::labs(x = "OR and 95% CI", y = "") +
    ggplot2::theme(panel.grid.major.y = ggplot2::element_blank(), 
                   panel.grid.minor = ggplot2::element_blank(),
                   axis.text = ggplot2::element_text(size=8),
                   text = ggplot2::element_text(size=8),
                   legend.position = "bottom",
                   legend.title = ggplot2::element_blank())
  
  ggplot2::ggsave(filename = paste0("output/mvmr_",i,".tiff"),
                  dpi = 300,width = 215.9, height = 279.4, unit = "cm", scale = 0.6)
  
}