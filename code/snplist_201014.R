rm(list=ls())
graphics.off()

# Specify paths ----------------------------------------------------------------

source("code/specify_paths.R")


# Load instruments -------------------------------------------------------------

instruments <- data.table::fread("data/instruments_all.txt", data.table = FALSE)

# Format exposure data ---------------------------------------------------------

exp <- TwoSampleMR::format_data(dat = instruments,
                                type = "exposure",
                                phenotype_col = "exposure")

# Clump data -------------------------------------------------------------------

exp <- TwoSampleMR::clump_data(exp,
                               clump_kb = 10000,
                               clump_r2 = 0.001,
                               pop = "EUR")

data.table::fwrite(exp,"data/instruments_clumped.txt")

# Format for sharing -----------------------------------------------------------

df <- exp[,c("SNP","chr.exposure","pos.exposure","effect_allele.exposure","other_allele.exposure")]
df <- unique(df)
colnames(df) <- c("SNP","chr","pos","effect_allele","other_allele")

data.table::fwrite(df,"output/snplist_201014.txt")

