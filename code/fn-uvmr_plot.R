uvmr_plot <- function(dat, type, trait) {
   
   # Record trait long name ----------------------------------------------------
   
   trait_long <- ifelse(trait=="cad","coronary artery disease",
                        ifelse(trait=="pad","peripheral artery disease",
                               ifelse(trait=="t2d","type 2 diabetes","")))

   # Restrict data to that needed for plotting ----------------------------------
   
   dat <- dat[,c("id.exposure","exposure","id.outcome","outcome","method","est","est_lci","est_uci","pval","nsnp","evidence")]
   dat <- unique(dat)
   
   # Convert data to wide (i.e. one row per exposure) ---------------------------
   
   dat$method <- ifelse(dat$method %in% c("Inverse variance weighted","Wald ratio"),"main",dat$method)
   
   dat_main <- dat[dat$method=="main",c("id.exposure","exposure","id.outcome","outcome","nsnp",
                                        "est","est_lci","est_uci","pval","evidence")]
   
   colnames(dat_main) <- c("id.exposure","exposure","id.outcome","outcome","nsnp",
                           "main_est","main_est_lci","main_est_uci","main_pval","evidence")
   
   dat_sens <- dat[dat$method!="main",c("id.exposure","exposure","id.outcome","outcome","nsnp",
                                        "method","est","est_lci","est_uci","pval")]
   
   colnames(dat_sens) <- c("id.exposure","exposure","id.outcome","outcome","nsnp",
                           "sensitivity_analysis","sensitivity_est","sensitivity_est_lci","sensitivity_est_uci","sensitivity_pval")
   
   dat <- merge(dat_main, dat_sens, by = c("id.exposure","exposure","id.outcome","outcome","nsnp"))
   
   if (type=="outcome") {
      
      # Format data ahead of plotting ----------------------------------------------
      
      dat$trait_plot <- paste0("OR: ",
                               sprintf("%.2f",dat$main_est),
                               "; 95% CI: ",
                               sprintf("%.2f",dat$main_est_lci),
                               " to ",
                               sprintf("%.2f",dat$main_est_uci),
                               "; ",
                               dat$nsnp," SNPs")
      
      # Specify breaks and labels --------------------------------------------------
      
      log_breaks <- 2^(seq(-5,15,1))
      log_labels <- c(sprintf("%.2f",2^(seq(-5,-2,1))),"0.5",sprintf("%.0f",2^(seq(0,15,1))))
      
      # Plot main results ----------------------------------------------------------
      
      main <- unique(dat[dat$evidence==TRUE,c("id.exposure","exposure","id.outcome","outcome","nsnp",
                            "main_est","main_est_lci","main_est_uci","main_pval",
                            "trait_plot")])
      
      main$yorder <- rank(main$main_est, ties.method = "random")
      
      ggplot2::ggplot(data = main,
                      mapping = ggplot2::aes(x = main_est, y = yorder)) +
         ggplot2::geom_vline(xintercept=1, col = "dark grey") +
         ggplot2::geom_linerange(ggplot2::aes(xmin = main_est_lci, xmax = main_est_uci), alpha = 0.5, size = 1, col = "#084594") +
         ggplot2::geom_point(shape = 15, size = 0.5, col = "#084594") +
         ggplot2::scale_x_continuous(trans = "log", 
                                     breaks = log_breaks,
                                     labels = log_labels,
                                     lim = c(round(min(main$main_est_lci),2)-0.01,round(max(main$main_est_uci),2)+1)) +
         ggplot2::scale_y_continuous(breaks = main$yorder,
                                     labels = main$exposure,
                                     sec.axis = ggplot2::sec_axis(~.,
                                                                  breaks = main$yorder,
                                                                  labels = main$trait_plot),
                                     expand = c(0,0.6)) +
         ggplot2::theme_minimal() +
         ggplot2::labs(x = paste0("OR and 95% confidence interval for\nthe effect of the risk factors on ",trait_long), 
                       y = "") +
         ggplot2::theme(panel.grid.major.y = ggplot2::element_blank(),
                        panel.grid.minor = ggplot2::element_blank(),
                        axis.text = ggplot2::element_text(size=8),
                        text = ggplot2::element_text(size=8))
      
      ggplot2::ggsave(filename = paste0("output/uvmr_",trait,"_",type,".jpeg"),
                      dpi = 300, width = 210, height = 297, unit = "mm", scale = 0.9)
      
      # Plot sensitivity analysis results ------------------------------------------
      
      ggplot2::ggplot(data = dat, mapping = ggplot2::aes(x = main_est, y = sensitivity_est)) +
         ggplot2::geom_hline(yintercept = 1, col = "dark gray") +
         ggplot2::geom_vline(xintercept = 1, col = "dark gray") +
         ggplot2::geom_abline(intercept = 0, slope = 1, col = "dark gray") +
         ggplot2::geom_errorbarh(ggplot2::aes(xmin = main_est_lci, xmax = main_est_uci), alpha = 0.5, col = "#084594") +
         ggplot2::geom_errorbar(ggplot2::aes(ymin = sensitivity_est_lci, ymax = sensitivity_est_uci), alpha = 0.5, col = "#084594") +
         ggplot2::scale_x_continuous(trans = "log",
                                     breaks = log_breaks,
                                     labels = log_labels,
                                     lim = c(0.03,32768)) +
         ggplot2::scale_y_continuous(trans = "log",
                                     breaks = log_breaks,
                                     labels = log_labels,
                                     lim = c(0.03,32768)) +
         ggplot2::theme_minimal() +
         ggplot2::labs(x = paste0("Estimate and 95% confidence interval for the effect of the risk factors on\n",trait_long," using the inverse variance weighted method"),
                       y = paste0("Estimate and 95% confidence interval for the effect of the risk factors on\n",trait_long," using the specified sensitivity analysis method")) +
         ggplot2::theme(panel.grid.major = ggplot2::element_blank(),
                        panel.grid.minor = ggplot2::element_blank(),
                        axis.text = ggplot2::element_text(size=8),
                        text = ggplot2::element_text(size=8),
                        legend.title = ggplot2::element_blank(),
                        legend.position = "bottom",
                        strip.placement = "outside",
                        axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust=1)) +
         ggplot2::facet_wrap(.~sensitivity_analysis, strip.position = "left", scales = "free")
      
      ggplot2::ggsave(filename = paste0("output/uvmr_",trait,"_",type,"_sa.jpeg"),
                      dpi = 300, height = 210, width = 210, unit = "mm", scale = 1)
      
   } else {
      
      # Format data ahead of plotting ----------------------------------------------
      
      dat$trait_plot <- paste0("Beta: ",
                               sprintf("%.2f",dat$main_est),
                               "; 95% CI: ",
                               sprintf("%.2f",dat$main_est_lci),
                               " to ",
                               sprintf("%.2f",dat$main_est_uci),
                               "; ",
                               dat$nsnp," SNPs")
      
      # Plot main results ----------------------------------------------------------
      
      main <- unique(dat[dat$evidence==TRUE,c("id.exposure","exposure","id.outcome","outcome","nsnp",
                            "main_est","main_est_lci","main_est_uci","main_pval",
                            "trait_plot")])
      
      main$yorder <- rank(main$main_est, ties.method = "random")
      
      ggplot2::ggplot(data = main, mapping = ggplot2::aes(x = main_est, y = yorder)) +
         ggplot2::geom_vline(xintercept=0, col = "dark grey") +
         ggplot2::geom_linerange(ggplot2::aes(xmin = main_est_lci, xmax = main_est_uci), alpha = 0.5, size = 1, col = "#084594") +
         ggplot2::geom_point(shape = 15, size = 0.5, col = "#084594") +
         ggplot2::scale_x_continuous(breaks = seq(-2,2,0.5),
                                     lim = c(-2,2)) +
         ggplot2::scale_y_continuous(breaks = main$yorder,
                                     labels = main$outcome,
                                     sec.axis = ggplot2::sec_axis(~.,
                                                                  breaks = main$yorder,
                                                                  labels = main$trait_plot),
                                     expand = c(0,0.6)) +
         ggplot2::theme_minimal() +
         ggplot2::labs(x = paste0("Beta and 95% confidence interval for\nthe effect of ",trait_long," on the risk factors"), 
                       y = "") +
         ggplot2::theme(panel.grid.major.y = ggplot2::element_blank(),
                        panel.grid.minor = ggplot2::element_blank(),
                        axis.text = ggplot2::element_text(size=8),
                        text = ggplot2::element_text(size=8))
      
      ggplot2::ggsave(filename = paste0("output/uvmr_",trait,"_",type,".jpeg"),
                      dpi = 300, width = 210, height = 297, unit = "mm", scale = 0.9)
      
      # Plot sensitivity analysis results ------------------------------------------
      
      ggplot2::ggplot(data = dat, mapping = ggplot2::aes(x = main_est, y = sensitivity_est)) +
         ggplot2::geom_hline(yintercept = 0, col = "dark gray") +
         ggplot2::geom_vline(xintercept = 0, col = "dark gray") +
         ggplot2::geom_abline(intercept = 0, slope = 1, col = "dark gray") +
         ggplot2::geom_errorbarh(ggplot2::aes(xmin = main_est_lci, xmax = main_est_uci), alpha = 0.5, col = "#084594") +
         ggplot2::geom_errorbar(ggplot2::aes(ymin = sensitivity_est_lci, ymax = sensitivity_est_uci), alpha = 0.5, col = "#084594") +
         ggplot2::scale_x_continuous(breaks = seq(-2,2,0.5),
                                     lim = c(-2,2)) +
         ggplot2::scale_y_continuous(breaks = seq(-2,2,0.5),
                                     lim = c(-2,2)) +
         ggplot2::theme_minimal() +
         ggplot2::labs(x = paste0("Estimate and 95% confidence interval for the effect of ",trait_long," on\nthe risk factors using the inverse variance weighted method"),
                       y = paste0("Estimate and 95% confidence interval for the effect of ",trait_long," on\nthe risk factors using the specified sensitivity analysis method")) +
         ggplot2::theme(panel.grid.major = ggplot2::element_blank(),
                        panel.grid.minor = ggplot2::element_blank(),
                        axis.text = ggplot2::element_text(size=8),
                        text = ggplot2::element_text(size=8),
                        legend.title = ggplot2::element_blank(),
                        legend.position = "bottom",
                        strip.placement = "outside",
                        axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust=1)) +
         ggplot2::facet_wrap(.~sensitivity_analysis, strip.position = "left", scales = "free")
      
      ggplot2::ggsave(filename = paste0("output/uvmr_",trait,"_",type,"_sa.jpeg"),
                      dpi = 300, height = 210, width = 210, unit = "mm", scale = 1)
      
   }
   
}