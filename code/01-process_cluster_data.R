# Collate basic cluster information
clusters <- data.frame(cluster_id = paste0("C",1:9),
                       weight_id = paste0("W",1:9),
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

# Save cluster info
data.table::fwrite(clusters,"data/cluster-info.csv")

# Set cluster threshold to 0.75 (as in Udler et al)
cluster_threshold <- 0.75

# Load variant data
variant_matrix <- readxl::read_excel("raw/L2EU.H.mat.9_T2D_71traits_259snps_samplesizefilter_v4.xlsx",
                                     range = "A1:K260")

# Collapse variant matrix
variant_list <- tidyr::gather(variant_matrix,weight_id,weight,W1:W9,factor_key = TRUE)

# Determine variants with weight > threshold
variant_list$meets_threshold <- ifelse(variant_list$weight>=cluster_threshold, TRUE, FALSE)

# Add cluster info
variant_list <- merge(variant_list,clusters,by=c("weight_id"),all.x = TRUE)

# Save variant list
data.table::fwrite(variant_list,"data/cluster-variants.csv")

# Load trait data
trait_matrix <- readxl::read_excel("raw/L2EU.H.mat.9_T2D_71traits_259snps_samplesizefilter_v4.xlsx",
                                   range = "M1:V143")
# Collapse trait matrix
trait_list <- tidyr::gather(trait_matrix,weight_id,weight,W1:W9,factor_key = TRUE)

# Determine traits with weight > threshold
trait_list$meets_threshold <- ifelse(trait_list$weight>=cluster_threshold, TRUE, FALSE)

# Add cluster info
trait_list <- merge(trait_list,clusters,by=c("weight_id"),all.x = TRUE)

# Save trait list
data.table::fwrite(trait_list,"data/cluster-traits.csv")