uvmr_plot <- function(dat, type, trait) {

   # Restrict data to that needed for plotting ----------------------------------

   dat <- dat[,c("exposure","exposure_long","outcome","outcome_long","method","or","lci_or","uci_or","pval","nsnp")]
   dat <- unique(dat)
   dat <- na.omit(dat)
   
   # Convert data to wide (i.e. one row per exposure) ---------------------------
   
   dat$method <- ifelse(dat$method %in% c("Inverse variance weighted","Wald ratio"),"main",dat$method)
   
   dat$method <- gsub(" ","_",dat$method)
   
   dat <- tidyr::pivot_wider(dat,
                                 id_cols = c("exposure","exposure_long","outcome","outcome_long","nsnp"),
                                 names_from = "method",
                                 names_sep = "-",
                                 values_from = c("or","lci_or","uci_or","pval"),
                                 values_fill = NA)
   
   # Convert data to long with main result recorded alongside sensitivity -------
   
   dat <- tidyr::pivot_longer(dat,
                                  cols = c(paste0("or-",c("MR_Egger","Simple_mode","Weighted_mode","Weighted_median")),
                                           paste0("lci_or-",c("MR_Egger","Simple_mode","Weighted_mode","Weighted_median")),
                                           paste0("uci_or-",c("MR_Egger","Simple_mode","Weighted_mode","Weighted_median")),
                                           paste0("pval-",c("MR_Egger","Simple_mode","Weighted_mode","Weighted_median"))),
                                  names_to = "sensitivity")
   
   # Convert or, lci and uci data back to wide ----------------------------------
   
   dat$stat <- gsub("-.*","",dat$sensitivity)
   
   dat$sensitivity <- gsub(".*-","",dat$sensitivity)
   
   dat <- tidyr::pivot_wider(dat,
                                 id_cols = c("exposure","exposure_long","outcome","outcome_long","nsnp","or-main","lci_or-main","uci_or-main","pval-main","sensitivity"),
                                 names_from = "stat", 
                                 values_from = "value")
   
   # Format data ahead of plotting ----------------------------------------------
   
   dat$sensitivity<- gsub("_"," ",dat$sensitivity)
   colnames(dat) <- c("exposure","exposure_long","outcome","outcome_long","nsnp","main_or","main_lci","main_uci","main_pval","sensitivity_analysis","sensitivity_or","sensitivity_lci","sensitivity_uci","sensitivity_pval")
   
   # dat$trait_plot <- paste0("OR: ",
   #                                 ifelse(nchar(signif(dat$main_or, digits = 3))<3, paste0(signif(dat$main_or, digits = 3),".0"),signif(dat$main_or, digits = 3)),
   #                                 "; 95% CI: ",
   #                                 ifelse(nchar(signif(dat$main_lci, digits = 3))<3, paste0(signif(dat$main_lci, digits = 3),".0"),signif(dat$main_lci, digits = 3)),
   #                                 " to ",
   #                                 ifelse(nchar(signif(dat$main_uci, digits = 3))<3, paste0(signif(dat$main_uci, digits = 3),".0"),signif(dat$main_uci, digits = 3)),
   #                                 "; ",
   #                                 dat$nsnp," SNPs")
   
   dat$trait_plot <- paste0("OR: ",
                            sprintf("%.2f",dat$main_or),
                            "; 95% CI: ",
                            sprintf("%.2f",dat$main_lci),
                            " to ",
                            sprintf("%.2f",dat$main_uci),
                            "; ",
                            dat$nsnp," SNPs")

    
   if (type=="exposure") {
     dat$x <- dat$outcome_long
   }
   
   if (type=="outcome") {
     dat$x <- dat$exposure_long
   }
  
   
   # Specify breaks and labels --------------------------------------------------
   
   log_breaks <- c(0.03,0.06,0.12,0.25,0.5,1,2,4,8,16,32,64,128,256)
   log_labels <- c("0.03","0.06","0.12","0.25","0.5","1","2","4","8","16","32","64","128","256")
   
   # Plot main results ----------------------------------------------------------
   
   main <- unique(dat[,c("exposure","exposure_long","outcome","outcome_long","nsnp",
                  "main_or","main_lci","main_uci","main_pval",
                  "trait_plot","x")])
  
   main$yorder <- rank(main$main_or, ties.method = "random")
   
ggplot2::ggplot(data = main,
                   mapping = ggplot2::aes(x = main_or, y = yorder)) +
     ggplot2::geom_vline(xintercept=1, col = "dark grey") +
     ggplot2::geom_linerange(ggplot2::aes(xmin = main_lci, xmax = main_uci), alpha = 0.5, size = 1, col = "#084594") +
     ggplot2::geom_point(shape = 15, size = 0.5, col = "#084594") +
     ggplot2::scale_x_continuous(trans = "log", 
                                 breaks = log_breaks,
                                 labels = log_labels,
                                 lim = c(round(min(main$main_lci),2)-0.01,round(max(main$main_uci),2)+1)) +
     ggplot2::scale_y_continuous(breaks = main$yorder,
                                 labels = main$x,
                                 sec.axis = ggplot2::sec_axis(~.,
                                                              breaks = main$yorder,
                                                              labels = main$trait_plot),
                                 expand = c(0,0.6)) +
     ggplot2::theme_minimal() +
     ggplot2::labs(x = ifelse(type=="outcome", 
                              paste0("Odds ratio and 95% confidence interval for\nthe effect of the features on ",main$outcome_long[1]),
                              paste0("Odds ratio and 95% confidence interval for\nthe effect of ",main$exposure_long[1]," on the features")), 
                   y = "") +
     ggplot2::theme(panel.grid.major.y = ggplot2::element_blank(),
                    panel.grid.minor = ggplot2::element_blank(),
                    axis.text = ggplot2::element_text(size=8),
                    text = ggplot2::element_text(size=8))
   
   ggplot2::ggsave(filename = paste0("output/uvmr_",trait,"_",type,".jpeg"),
                   dpi = 300, width = 210, height = 297, unit = "mm", scale = 0.9)
   
   # Plot sensitivity analysis results ------------------------------------------
   
   ggplot2::ggplot(data = dat[dat$nsnp>9,],
                   mapping = ggplot2::aes(x = main_or, y = sensitivity_or)) +
     ggplot2::geom_hline(yintercept = 1, col = "dark gray") +
     ggplot2::geom_vline(xintercept = 1, col = "dark gray") +
     ggplot2::geom_abline(intercept = 0, slope = 1, col = "dark gray") +
     ggplot2::geom_errorbarh(ggplot2::aes(xmin = main_lci, xmax = main_uci), alpha = 0.5, col = "#084594") +
     ggplot2::geom_errorbar(ggplot2::aes(ymin = sensitivity_lci, ymax = sensitivity_uci), alpha = 0.5, col = "#084594") +
     ggplot2::scale_x_continuous(trans = "log", 
                                 breaks = log_breaks,
                                 labels = log_labels,
                                 lim = c(round(min(dat[dat$nsnp>9,]$main_lci,dat[dat$nsnp>9,]$sensitivity_lci),2)-0.01,round(max(dat[dat$nsnp>9,]$main_uci,dat[dat$nsnp>9,]$sensitivity_uci),2)+1)) +
     ggplot2::scale_y_continuous(trans = "log", 
                                 breaks = log_breaks,
                                 labels = log_labels,
                                 lim = c(round(min(dat[dat$nsnp>9,]$main_lci,dat[dat$nsnp>9,]$sensitivity_lci),2)-0.01,round(max(dat[dat$nsnp>9,]$main_uci,dat[dat$nsnp>9,]$sensitivity_uci),2)+1)) +
     ggplot2::theme_minimal() +
     ggplot2::labs(x = ifelse(type=="outcome", 
                              paste0("Odds ratio and 95% confidence interval for the effect of the features on\n",dat$outcome_long[1]," using the inverse variance weighted method"),
                              paste0("Odds ratio and 95% confidence interval for the effect of ",dat$exposure_long[1]," on\nthe features using the inverse variance weighted method")),
                   y = ifelse(type=="outcome",
                             paste0("Odds ratio and 95% confidence interval for the effect of the features on\n",dat$outcome_long[1]," using the specified sensitivity analysis method"),
                             paste0("Odds ratio and 95% confidence interval for the effect of ",dat$exposure_long[1]," on\nthe features using the specified sensitivity analysis method"))) +
     ggplot2::theme(panel.grid.major = ggplot2::element_blank(),
                    panel.grid.minor = ggplot2::element_blank(),
                    axis.text = ggplot2::element_text(size=8),
                    text = ggplot2::element_text(size=8),
                    legend.title = ggplot2::element_blank(),
                    legend.position = "bottom",
                    strip.placement = "outside") +
     ggplot2::facet_wrap(.~sensitivity_analysis, strip.position = "left", scales = "free")
   
   ggplot2::ggsave(filename = paste0("output/uvmr_",trait,"_",type,"_sa.jpeg"),
                   dpi = 300, height = 210, width = 210, unit = "mm", scale = 1)
   
 }