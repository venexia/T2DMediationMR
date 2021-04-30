rm(list=ls())
graphics.off()

# Load UVMR results ------------------------------------------------------------

uvmr <- data.table::fread("output/results.csv",
                        stringsAsFactors = FALSE,
                        data.table = FALSE)

# Add confidence intervals to UVMR results -------------------------------------

uvmr$lci <- uvmr$b - qnorm(0.975)*uvmr$se
uvmr$uci <- uvmr$b + qnorm(0.975)*uvmr$se

# Convert UVMR results to odds ratios ------------------------------------------

uvmr$or <- exp(uvmr$b)
uvmr$lci_or <- exp(uvmr$lci)
uvmr$uci_or <- exp(uvmr$uci)

# Restrict UVMR results to relevant results ------------------------------------ 

uvmr <- uvmr[uvmr$method %in% c("Wald ratio","Inverse variance weighted"),]
uvmr <- uvmr[uvmr$exposure %in% c("t2d","t2d_linear") | uvmr$outcome %in% c("t2d","t2d_linear"),]
uvmr <- uvmr[,c("exposure", "outcome", "or", "lci_or", "uci_or")]

# Identify linear and logistic UVMR results ------------------------------------

uvmr$model <- ifelse(grepl("t2d",uvmr$exposure),uvmr$exposure,uvmr$outcome)
uvmr$model <- ifelse(uvmr$model=="t2d_linear","linear","logistic")

uvmr$exposure <- ifelse(uvmr$exposure=="t2d_linear","t2d",uvmr$exposure)
uvmr$outcome <- ifelse(uvmr$outcome=="t2d_linear","t2d",uvmr$outcome)

uvmr <- tidyr::pivot_wider(uvmr, names_from = "model", values_from = c("or","lci_or","uci_or"))
uvmr$effect <- "total"

# Load MVMR results ------------------------------------------------------------

mvmr <- data.table::fread("output/twostep_results.csv", data.table = FALSE)
mvmr$model <- ifelse(mvmr$mediator=="t2d_linear","linear","logistic")

# Add confidence intervals to MVMR results -------------------------------------

mvmr$lci <- mvmr$estimate - qnorm(0.975)*mvmr$se
mvmr$uci <- mvmr$estimate + qnorm(0.975)*mvmr$se

# Convert MVMR results to odds ratios ------------------------------------------

mvmr$or <- exp(mvmr$estimate)
mvmr$lci_or <- exp(mvmr$lci)
mvmr$uci_or <- exp(mvmr$uci)

# Identify linear and logistic MVMR results ------------------------------------

mvmr <- mvmr[mvmr$effect %in% c("direct","indirect"), c("exposure","outcome","effect","model","or","lci_or","uci_or")]
mvmr <- tidyr::pivot_wider(mvmr, names_from = "model", values_from = c("or","lci_or","uci_or"))

# Combine UVMR and MVMR results ------------------------------------------------

df <- rbind(mvmr,uvmr)

# Plot linear vs logistic comparison -------------------------------------------

ggplot2::ggplot(df, mapping = ggplot2::aes(x = or_linear, y = or_logistic, color = effect)) +
  ggplot2::geom_hline(yintercept = 1, col = "dark gray") +
  ggplot2::geom_vline(xintercept = 1, col = "dark gray") +
  ggplot2::geom_abline(intercept = 0, slope = 1, col = "dark gray") +
  ggplot2::geom_point(alpha = 0.5) +
  ggplot2::scale_x_continuous(trans = "log", 
                              breaks = c(0.5,1,2,4,8),
                              lim = c(0.5,9)) +
  ggplot2::scale_y_continuous(trans = "log", 
                              breaks = c(0.5,1,2,4,8),
                              lim = c(0.5,9)) +
  ggplot2::theme_minimal() +
  ggplot2::labs(x = "Mendelian randomization estimate using the linear type 2 diabetes GWAS", 
                y = "Mendelian randomization estimate using the logistic type 2 diabetes GWAS") +
  ggplot2::scale_color_manual(breaks = c("total","direct","indirect"), 
                                values = c("#377eb8","#984ea3","#e41a1c"), 
                                labels = c("Univariate MR","Two-step MR for mediation, direct effect","Two-step MR for mediation, indirect effect"),
                                name = "") +
  ggplot2::theme(panel.grid.major = ggplot2::element_blank(),
                 panel.grid.minor = ggplot2::element_blank(),
                 axis.text = ggplot2::element_text(size=8),
                 text = ggplot2::element_text(size=8),
                 legend.position = "bottom")

# Save linear vs logistic comparison plot --------------------------------------

ggplot2::ggsave(filename = paste0("output/linear_logistic_comparison.jpeg"),
                dpi = 300, height = 210, width = 210, unit = "mm", scale = 1)
