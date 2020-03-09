# Note basic cluster information

clusters <- data.frame(cluster_id = paste0("C",1:9),
                       cluster_name = c("Proinsulin",
                                        "Lipodystrophy",
                                        "LipoproteinA",
                                        "Liver_lipid",
                                        "Beta-cell2",
                                        "Beta-cell",
                                        "Obesity",
                                        "Blood_traits",
                                        "Cholesterol"),
                       stringsAsFactors = FALSE)

for (i in 1:nrow(clusters)) {
  
  # Load variant data
  tmp <- read_excel("raw/L2EU.H.mat.9_T2D_71traits_259snps_samplesizefilter_v4.xlsx",
                    sheet = clusters$cluster_name[i], range = "A1:K260")
  
  # Save variant data
  data.table::fwrite(tmp,paste0("data/",clusters$cluster_id[i],"-variants"))
  
  # Load trait data
  tmp <- read_excel("raw/L2EU.H.mat.9_T2D_71traits_259snps_samplesizefilter_v4.xlsx",
                    sheet = clusters$cluster_name[i], range = "M1:V143")
  
  # Save trait data
  data.table::fwrite(tmp,paste0("data/",clusters$cluster_id[i],"-traits"))
  
}