rm(list=ls())
graphics.off()

# Load packages ----------------------------------------------------------------

library(tidyverse)

# Specify paths ----------------------------------------------------------------

source("code/specify_paths.R")

# Source functions -------------------------------------------------------------

source("code/fn-rsid.R")

# Keep list --------------------------------------------------------------------

keep <- c("keep","clean_gwas","MVMR","rsid",ls(pattern="path_"))

