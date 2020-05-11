rm(list=ls())
graphics.off()

# Load packages ----------------------------------------------------------------

library(tidyverse)

# Specify paths ----------------------------------------------------------------

source("code/specify_paths.R")

# Source functions -------------------------------------------------------------

source("code/fn-rsid.R")
source("code/fn-extract_feature_ieugwas.R")
source("code/fn-extract_feature_manual.R")
source("code/fn-create_analysis_data.R")

# Keep list --------------------------------------------------------------------

keep <- c("keep","clean_gwas","MVMR","rsid","create_analysis_data","extract_feature_ieugwas","extract_feature_manual",ls(pattern="path_"))

