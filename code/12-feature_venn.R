rm(list=ls())
graphics.off()

# Specify paths ----------------------------------------------------------------

source("code/specify_paths.R", echo = TRUE)

# Load evidence ----------------------------------------------------------------

features <- data.table::fread("output/evidence_summary.csv", data.table = FALSE)

# Annonate features ------------------------------------------------------------

gwas <- data.table::fread("raw/gwas.csv",
                          stringsAsFactors = FALSE,
                          select = c("trait","trait_long"),
                          data.table = FALSE)

features <- merge(features, gwas, by.x = "feature", by.y = "trait", all.x = TRUE)

# Record features not in Venn diagram ------------------------------------------

noVenn <- features[features$feature_t2d_evidence==FALSE &
                     features$feature_pad_evidence==FALSE &
                     features$feature_cad_evidence==FALSE &
                     features$t2d_feature_evidence==FALSE,]

paste0(nrow(noVenn)," of ",nrow(features),
       " features did not pass the FDR threshold of 0.05 for any analysis: ",
       paste(noVenn$trait_long, collapse = ", "))

# Generate Venn diagram --------------------------------------------------------

x <- list(feature_t2d_evidence = features[features$feature_t2d_evidence==TRUE,]$feature,
          feature_pad_evidence = features[features$feature_pad_evidence==TRUE,]$feature,
          feature_cad_evidence = features[features$feature_cad_evidence==TRUE,]$feature,
          t2d_feature_evidence = features[features$t2d_feature_evidence==TRUE,]$feature)

ggVennDiagram::ggVennDiagram(x, 
                             category.names = c("Feature > T2D" , "Feature > PAD", "Feature > CAD", "T2D > Feature"),
                             label = "count") +
  ggplot2::scale_fill_gradient(low = "white", high = "white") +
  ggplot2::theme(legend.position = "none")
  
  ggplot2::ggsave(filename = "output/feature_venn.jpeg",
                  dpi = 300, height = 210, width = 210, unit = "mm", scale = 1)

  ggVennDiagram::ggVennDiagram(x, 
                               category.names = c("Feature > T2D" , "Feature > PAD", "Feature > CAD", "T2D > Feature"),
                               label = NULL) +
    ggplot2::scale_fill_gradient(low = "white", high = "white") +
    ggplot2::theme(legend.position = "none")
  
  ggplot2::ggsave(filename = "output/feature_venn_nolab.jpeg",
                  dpi = 300, height = 210, width = 210, unit = "mm", scale = 1)
  