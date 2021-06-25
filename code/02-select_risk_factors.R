rm(list=ls())
graphics.off()

# Load packages ----------------------------------------------------------------

library(magrittr)

# Specify paths ----------------------------------------------------------------

source("code/specify_paths.R", echo = TRUE)

# List available outcomes ------------------------------------------------------

ao <- TwoSampleMR::available_outcomes()

# Select relevant variables ----------------------------------------------------

df <- ao[,c("id","trait","population","sex","consortium","author","category","ncase","sample_size","nsnp","note")]

# Remove proteins, expression, metabolomics, non-Neale UKB and raw GWAS --------

df <- df[!grepl("prot",df$id),]
df <- df[!grepl("eqtl",df$id),]
df <- df[!grepl("met",df$id),]
df <- df[!grepl("ukb-a",df$id),]
df <- df[!grepl("ukb-b",df$id),]
df <- df[!grepl("ukb-e",df$id),]
df <- df[!grepl("raw",df$id),]

# Remove adjusted GWAS ---------------------------------------------------------

df <- df[!grepl("adjust",tolower(df$note)),]
df <- df[df$consortium!="International Consortium of Blood Pressure",] # BMI adjusted
df <- df[df$consortium!="HaemGen",] # Prefer UKB

# Select correct UKB GWAS ------------------------------------------------------

df <- df[df$consortium!="UK Biobank",] # Non-Neale UKB GWAS
df <- rbind(df,ao[ao$id %in% c("ukb-b-10787","ukb-b-20175","ukb-b-7992","ukb-b-8909","ukb-b-15590","ukb-b-9405"),c("id","trait","population","sex","consortium","author","category","ncase","sample_size","nsnp","note")])
df$sample_size <- ifelse(is.na(df$sample_size) & grepl("ukb-d",df$id),361194,df$sample_size)

# Restricted based on population, sex, category, continous/binary --------------

df <- df[df$population %in% c("European","Mixed"),]
df <- df[!(df$sex %in% c("Males","Females")),]
df <- df[df$category %in% c("Continuous","Risk factor","Metabolites"),]
df <- df[is.na(df$ncase),]

# Remove irrelevant traits -----------------------------------------------------

df$trait <- tolower(df$trait)
df <- df[!grepl("age ",df$trait),]
df <- df[!grepl("years of ",df$trait),]
df <- df[!grepl("difference in ",df$trait),]
df <- df[!grepl("digit symbol",df$trait),]
df <- df[!grepl("home location",df$trait),]
df <- df[!grepl("road",df$trait),]
df <- df[!grepl("work",df$trait),]
df <- df[!grepl("inspection time",df$trait),]
df <- df[!grepl("reaction time",df$trait),]
df <- df[!grepl("symbol search",df$trait),]
df <- df[!grepl("night shifts",df$trait),]

# Correct trait names ----------------------------------------------------------

df$trait <- ifelse(df$trait=="cigarettes smoked per day", "cigarettes per day", df$trait)
df$trait <- ifelse(df$trait=="c-reactive protein level", "c-reactive protein", df$trait)
df$trait <- ifelse(df$trait=="glycated haemoglobin", "hba1c", df$trait)
df$trait <- ifelse(df$trait=="ldl direct", "ldl cholesterol", df$trait)
df$trait <- ifelse(df$trait=="serum creatinine (egfrcrea)", "creatinine", df$trait)
df$trait <- ifelse(df$trait=="transferrin saturation", "transferrin", df$trait)

# Retain largest GWAS for trait ------------------------------------------------

df <- df %>%
  dplyr::group_by(trait) %>%
  dplyr::top_n(n = 1, wt = sample_size)

# For GWAS of same size, retain GWAS with most SNPs ----------------------------

df <- df %>%
  dplyr::group_by(trait) %>%
  dplyr::top_n(n = 1, wt = nsnp)

# Format risk factor data for analysis -----------------------------------------

risk_factors <- ao[ao$id %in% df$id,c("trait","id","year","author","consortium","sex","pmid","population","sample_size","nsnp","build","unit","ncase","ncontrol")]
risk_factors$sample_size <- ifelse(is.na(risk_factors$sample_size) & grepl("ukb-d",risk_factors$id),361194,risk_factors$sample_size)
risk_factors$consortium <- ifelse(risk_factors$author=="Neale lab", "UK Biobank", risk_factors$consortium)
risk_factors$consortium <- ifelse(risk_factors$consortium=="NA", "", risk_factors$consortium)
risk_factors$type <- "risk_factor"

# Save risk factor data for analysis -------------------------------------------

data.table::fwrite(risk_factors,"data/risk_factors.csv")
