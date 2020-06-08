rm(list=ls())
graphics.off()

# Specify paths ----------------------------------------------------------------

source("code/specify_paths.R")

# Load feature data ------------------------------------------------------------

gwas <- data.table::fread("raw/gwas.csv",
                          stringsAsFactors = FALSE,
                          select = c("trait","trait_long"),
                          data.table = FALSE)

colnames(gwas) <- c("exposure","exposure_long")

features <- gwas[!(gwas$exposure %in% c("t2d","pad","cad")),]$exposure

# Load UVMR results -------------------------------------------------------------

df <- data.table::fread("output/uvmr_results.csv",
                        stringsAsFactors = FALSE,
                        data.table = FALSE)

# Annotate results -------------------------------------------------------------

df <- merge(gwas, df, by = c("exposure"))

# Add confidence intervals -----------------------------------------------------

df$lci <- df$b - qnorm(0.975)*df$se
df$uci <- df$b + qnorm(0.975)*df$se

# Convert to odds ratios -------------------------------------------------------

df$or <- exp(df$b)
df$lci_or <- exp(df$lci)
df$uci_or <- exp(df$uci)

# Make all odds ratios > 1 -----------------------------------------------------

df$or_pos <- ifelse(sign(df$b)==1,df$or,1/abs(df$or))
df$lci_or_pos <- ifelse(sign(df$b)==1,df$lci_or,1/abs(df$lci_or))
df$uci_or_pos <- ifelse(sign(df$b)==1,df$uci_or,1/abs(df$uci_or))

# Record direction of effect ---------------------------------------------------

df$direction <- ""
df$direction <- ifelse(sign(df$b)==1,"Increases risk",df$direction)
df$direction <- ifelse(sign(df$b)==-1,"Decreases risk",df$direction)

# Mark estimates that exclude the null -----------------------------------------

df$excl_null <- FALSE
df$excl_null <- ifelse(df$lci_or_pos>1 & df$uci_or_pos>1,TRUE,df$excl_null)
df$excl_null <- ifelse(df$lci_or_pos<1 & df$uci_or_pos<1,TRUE,df$excl_null)

# Rank estimates by size -------------------------------------------------------

tmp <- df[df$outcome=="t2d" & df$method %in% c("Wald ratio","Inverse variance weighted") & df$exposure %in% features,c("exposure","or_pos")]
tmp$order <- rank(-tmp$or_pos,na.last = TRUE)
tmp$or_pos <- NULL
data.table::fwrite(tmp,"output/t2d_assoc_rank.csv",row.names = FALSE)
df <- merge(df, tmp, by = c("exposure"), all.x = TRUE)

# Trim exposure labels ---------------------------------------------------------

df$exposure_long <- gsub("adjusted","adj.",df$exposure_long)
df$exposure_long <- gsub("distribution","dist.",df$exposure_long)
df$exposure_long <- gsub("concentration","con.",df$exposure_long)
df$exposure_long <- gsub("females","F",df$exposure_long)
df$exposure_long <- gsub("males","M",df$exposure_long)
df$exposure_long <- factor(df$exposure_long)

# Plot -------------------------------------------------------------------------

for (i in c("t2d","cad","pad")) {
  
  ggplot2::ggplot(data = df[df$outcome==i & df$method %in% c("Wald ratio","Inverse variance weighted") & df$exposure %in% features,], 
                  mapping = ggplot2::aes(x = forcats::fct_reorder(exposure_long, order, .desc = FALSE), y = or_pos, color = direction)) +
    ggplot2::geom_hline(yintercept=1, col = "dark gray") +
    ggplot2::geom_point(shape = 15, size = 0.5) +
    ggplot2::geom_linerange(ggplot2::aes(ymin = lci_or_pos, ymax = uci_or_pos), alpha = 0.6, size = 1) +
    ggplot2::scale_y_continuous(trans = "log", breaks = c(0.12,0.25,0.5,1,2,4,8,16),lim = c(0.06,32)) +
    ggplot2::scale_color_manual(breaks = c("Decreases risk","Increases risk"), values=c("#377eb8","#e41a1c")) +
    ggplot2::theme_minimal() +
    ggplot2::labs(y = "OR and 95% CI", x = "") +
    ggplot2::theme(panel.grid.major.x = ggplot2::element_blank(), 
                   panel.grid.minor = ggplot2::element_blank(),
                   axis.text.x = ggplot2::element_text(angle = 90,hjust = 1,vjust = 0.5),
                   axis.text = ggplot2::element_text(size=8),
                   text = ggplot2::element_text(size=8),
                   legend.position = "bottom",
                   legend.title = ggplot2::element_blank())
  
  ggplot2::ggsave(filename = paste0("output/uvmr_",i,".tiff"),
                  dpi = 300,width = 33.87, height = 19.05, unit = "cm", scale = 0.6)
  
}