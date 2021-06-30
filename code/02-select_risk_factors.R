rm(list=ls())
graphics.off()

# Load packages ----------------------------------------------------------------

library(magrittr)

# Specify paths ----------------------------------------------------------------

source("code/specify_paths.R", echo = TRUE)

# List available outcomes ------------------------------------------------------

ao <- TwoSampleMR::available_outcomes()

# Record number of traits ------------------------------------------------------

N <- data.frame(N = numeric(),
                criteria = character())

# Select relevant variables ----------------------------------------------------

df <- ao[,c("id","trait","population","sex","consortium","author","category","ncase","sample_size","nsnp","note")]
N[nrow(N)+1,] <- c(nrow(df),"IEU OpenGWAS database")

# Remove proteins, expression, metabolomics, non-Neale UKB and raw GWAS --------

df <- df[!grepl("eqtl",df$id),]
N[nrow(N)+1,] <- c(nrow(df),"Remove eQTLs")

df <- df[!grepl("prot",df$id),]
N[nrow(N)+1,] <- c(nrow(df),"Remove pQTLs")

df <- df[!grepl("met",df$id),]
N[nrow(N)+1,] <- c(nrow(df),"Remove mQTLs")

df <- df[!grepl("ukb-b",df$id),]
df <- df[!grepl("ukb-e",df$id),]
df <- df[!grepl("raw",df$id),]
df <- df[df$consortium!="UK Biobank",] # Non-Neale UKB GWAS
N[nrow(N)+1,] <- c(nrow(df),"Remove UK Biobank GWAS not from Neale lab")

# Remove adjusted GWAS ---------------------------------------------------------

df <- df[!grepl("adjust",tolower(df$note)),]
df <- df[df$consortium!="International Consortium of Blood Pressure",] 
df <- df[df$consortium!="HaemGen",]
N[nrow(N)+1,] <- c(nrow(df),"Remove adjusted GWAS")

# Correct UKB GWAS info --------------------------------------------------------

df$sample_size <- ifelse(is.na(df$sample_size) & grepl("ukb-d",df$id),361194,df$sample_size)
df$category <- ifelse(df$category=="NA" & grepl("ukb-a",df$id),"Continuous",df$category)

# Restricted based on population, sex, category, continous/binary --------------

df <- df[df$population %in% c("European","Mixed"),]
N[nrow(N)+1,] <- c(nrow(df),"Restrict to European or mixed ancestry GWAS")

df <- df[!(df$sex %in% c("Males","Females")),]
N[nrow(N)+1,] <- c(nrow(df),"Remove male only and female only GWAS")

df <- df[df$category %in% c("Continuous","Risk factor","Metabolites"),]
df <- df[is.na(df$ncase),]
N[nrow(N)+1,] <- c(nrow(df),"Remove GWAS of binary traits")

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
df <- df[!grepl("time spent",df$trait),]
df <- df[!grepl("frequency of",df$trait),]
df <- df[!grepl("number of",df$trait),]
df <- df[!grepl("preview only",df$trait),]
N[nrow(N)+1,] <- c(nrow(df),"Remove non-biological traits")

df <- df[!grepl("(left)",df$trait),]
df <- df[!grepl("(right)",df$trait),]
N[nrow(N)+1,] <- c(nrow(df),"Remove left and right specific traits")

# Correct trait names ----------------------------------------------------------

df$trait <- ifelse(df$trait=="standing height", "height", df$trait)
df$trait <- ifelse(df$trait=="cigarettes smoked per day", "cigarettes per day", df$trait)
df$trait <- ifelse(df$trait=="c-reactive protein level", "c-reactive protein", df$trait)
df$trait <- ifelse(df$trait=="glycated haemoglobin", "hba1c", df$trait)
df$trait <- ifelse(df$trait=="ldl direct", "ldl cholesterol", df$trait)
df$trait <- ifelse(df$trait=="serum creatinine (egfrcrea)", "creatinine", df$trait)
df$trait <- ifelse(df$trait=="transferrin saturation", "transferrin", df$trait)
df$trait <- ifelse(df$trait=="weight, manual entry", "weight", df$trait)
df$trait <- ifelse(df$trait=="body mass index (bmi)", "body mass index", df$trait)
df$trait <- ifelse(df$trait=="Creatinine (enzymatic) in urine", "Creatinine", df$trait)
df$trait <- gsub("  best measure","",df$trait)
df$trait <- gsub("  predicted percentage","",df$trait)
df$trait <- gsub("  predicted","",df$trait)
df$trait <- gsub("  automated reading","",df$trait)

# Retain largest GWAS for trait ------------------------------------------------

df <- df %>%
  dplyr::group_by(trait) %>%
  dplyr::top_n(n = 1, wt = sample_size)

N[nrow(N)+1,] <- c(nrow(df),"Retain largest GWAS for trait by sample size")

# For GWAS of same size, retain GWAS with most SNPs ----------------------------

df <- df %>%
  dplyr::group_by(trait) %>%
  dplyr::top_n(n = 1, wt = nsnp)

N[nrow(N)+1,] <- c(nrow(df),"Retain GWAS with most SNPs for trait when sample size is equal")

# Format risk factor data for analysis -----------------------------------------

risk_factors <- ao[ao$id %in% df$id,c("trait","id","year","author","consortium","sex","pmid","population","sample_size","nsnp","build","unit","ncase","ncontrol")]
risk_factors$sample_size <- ifelse(is.na(risk_factors$sample_size) & grepl("ukb-d",risk_factors$id),361194,risk_factors$sample_size)
risk_factors$consortium <- ifelse(risk_factors$author=="Neale lab", "UK Biobank", risk_factors$consortium)
risk_factors$consortium <- ifelse(risk_factors$consortium=="NA", "", risk_factors$consortium)
risk_factors$type <- "risk_factor"

# Save risk factor data for analysis -------------------------------------------

data.table::fwrite(risk_factors,"data/risk_factors.csv")
data.table::fwrite(N, "output/risk_factor_selection.csv")
