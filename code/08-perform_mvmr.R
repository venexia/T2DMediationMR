rm(list=ls())
graphics.off()

# Specify paths ----------------------------------------------------------------

source("code/specify_paths.R", echo = TRUE)

# Source functions -------------------------------------------------------------

source("code/fn-mvmr.R", echo = TRUE)

# Load instruments -------------------------------------------------------------

instruments <- data.table::fread("data/instruments.csv", data.table = FALSE)

# Load evidence summary --------------------------------------------------------

evidence <- data.table::fread("output/evidence_summary.csv", data.table = FALSE)
evidence <- evidence[(evidence$rf_t2d_evidence+evidence$t2d_rf_evidence)==1 & (evidence$rf_pad_evidence+evidence$rf_cad_evidence)>0,]
evidence <- evidence[!is.na(evidence$trait),]

# Determine analyses to be run -------------------------------------------------

analyses <- evidence[,c("id","rf_t2d_evidence","t2d_rf_evidence","rf_pad_evidence","rf_cad_evidence")]
analyses <- tidyr::pivot_longer(analyses, cols = c("rf_pad_evidence","rf_cad_evidence"), names_to = "outcome")
analyses <- analyses[analyses$value==TRUE, c("id","rf_t2d_evidence","t2d_rf_evidence","outcome")]
analyses$outcome <- gsub("rf_","",gsub("_evidence","",analyses$outcome))
analyses$exposure <- ifelse(analyses$rf_t2d_evidence==TRUE & analyses$t2d_rf_evidence==FALSE, analyses$id, "t2d")
analyses$mediator <- ifelse(analyses$rf_t2d_evidence==FALSE & analyses$t2d_rf_evidence==TRUE, analyses$id, "t2d")
analyses <- analyses[,c("exposure","mediator","outcome")]
tmp <- analyses
tmp$exposure <- ifelse(analyses$exposure=="t2d","t2d_linear",analyses$exposure)
tmp$mediator <- ifelse(analyses$mediator=="t2d","t2d_linear",analyses$mediator)
analyses <- rbind(analyses,tmp)
analyses$analysis <- paste(analyses$exposure, analyses$mediator, analyses$outcome, sep = "/")
data.table::fwrite(analyses,"data/mvmr_summary.csv")

# Load source data info --------------------------------------------------------

risk_factors <- data.table::fread("data/risk_factors.csv", data.table = FALSE)

# Perform analyses -------------------------------------------------------------

results <- data.table::fread("output/mvmr_results.csv")

for (i in 29:nrow(analyses)) {
  
  tmp_results <- mvmr(exposure = analyses$exposure[i],
                      mediator = analyses$mediator[i],
                      outcome = analyses$outcome[i])
  
  results <- rbind(results, tmp_results)
  
}

# Save results -----------------------------------------------------------------

data.table::fwrite(results,"output/mvmr_results.csv",row.names = FALSE)