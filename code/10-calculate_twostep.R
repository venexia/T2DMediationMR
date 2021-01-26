rm(list=ls())
graphics.off()

# Load filtered features list --------------------------------------------------

features <- data.table::fread("output/evidence_summary.csv", data.table = FALSE)

features$exposure <- NA
features$exposure <- ifelse(features$feature_t2d_evidence==TRUE & features$t2d_feature_evidence==FALSE,features$feature,features$exposure)
features$exposure <- ifelse(features$feature_t2d_evidence==FALSE & features$t2d_feature_evidence==TRUE,"t2d",features$exposure)

features$mediator <- NA
features$mediator <- ifelse(features$feature_t2d_evidence==TRUE & features$t2d_feature_evidence==FALSE,"t2d",features$mediator)
features$mediator <- ifelse(features$feature_t2d_evidence==FALSE & features$t2d_feature_evidence==TRUE,features$feature,features$mediator)

features <- tidyr::pivot_longer(features,
                                cols = c("feature_pad_evidence","feature_cad_evidence"),
                                names_to = "analysis",
                                values_to = "evidence")

features$outcome <- NA
features$outcome <- ifelse(features$analysis=="feature_pad_evidence" & features$evidence==TRUE,"pad",features$outcome)
features$outcome <- ifelse(features$analysis=="feature_cad_evidence" & features$evidence==TRUE,"cad",features$outcome)

features <- na.omit(features)
features <- features[,c("exposure","mediator","outcome")]

# Prepare MVMR results ---------------------------------------------------------

mvmr <-  data.table::fread("output/mvmr_results.csv",
                           select = c("analysis","exposure","outcome","effect","estimate","se"),
                           stringsAsFactors = FALSE,
                           data.table = FALSE)

mvmr <- mvmr[mvmr$effect=="direct",]
mvmr$effect <- NULL
colnames(mvmr) <- c("analysis","mediator","outcome","estimate_B","se_B")

mvmr$exposure <- substr(mvmr$analysis,1,nchar(mvmr$analysis)-4)
mvmr$exposure <- ifelse(mvmr$exposure==mvmr$mediator,"t2d",mvmr$exposure)
mvmr$analysis <- NULL

mvmr <- merge(features,mvmr,all.x = TRUE, by = c("exposure","mediator","outcome"))
mvmr$exposure_mediator <- paste0(mvmr$exposure,"_",mvmr$mediator)

# Prepare UVMR results ---------------------------------------------------------

uvmr <-  data.table::fread("output/results.csv",
                           select = c("exposure","outcome","method","b","se"),
                           stringsAsFactors = FALSE,
                           data.table = FALSE)

uvmr <- uvmr[uvmr$method %in% c("Inverse variance weighted","Wald ratio"), c("exposure","outcome","b","se")]

colnames(uvmr) <- c("exposure","mediator","estimate_A","se_A")
uvmr$exposure_mediator <- paste0(uvmr$exposure,"_",uvmr$mediator)

uvmr <- uvmr[uvmr$exposure_mediator %in% mvmr$exposure_mediator,]

# Combine estimates ------------------------------------------------------------

df <- merge(uvmr, mvmr, by = c("exposure","mediator","exposure_mediator"), all = TRUE)
df$exposure_mediator <- NULL
df <- na.omit(df)

# Calculate indirect effects ---------------------------------------------------

df$estimate <- df$estimate_A * df$estimate_B
# df$se <- abs(df$estimate) * sqrt((df$se_A/df$estimate_A)^2 + (df$se_B/df$estimate_B)^2)
df$se <- sqrt((df$se_A)^2 + (df$se_B)^2)

df <- df[,c("exposure","mediator","outcome","estimate_A","se_A","estimate_B","se_B","estimate","se")]
df$pc_A <- df$se_A/abs(df$estimate_A)
df$pc_B <- df$se_B/abs(df$estimate_B)
df$pc <- df$se/abs(df$estimate)

df <- df[,c("exposure","mediator","outcome","estimate","se")]
df$effect <- "indirect"
df$analysis <- paste0(df$exposure,"_",df$mediator,"_",df$outcome)

# Add direct effect results ----------------------------------------------------

direct <-  data.table::fread("output/mvmr_results.csv",
                             select = c("analysis","exposure","outcome","effect","estimate","se"),
                             stringsAsFactors = FALSE,
                             data.table = FALSE)

direct$mediator <- ifelse(direct$exposure=="t2d", gsub("_.*","",direct$analysis), "t2d")
direct$analysis <- paste0(direct$exposure,"_",direct$mediator,"_",direct$outcome)

direct <- direct[direct$analysis %in% df$analysis & direct$effect=="direct", 
                 c("exposure","mediator","outcome","estimate","se","effect","analysis")]

df <- rbind(df,direct)

# Add total effect results -----------------------------------------------------

df$exp_out <- paste0(df$exposure,"_",df$outcome)

total <-  data.table::fread("output/results.csv",
                            stringsAsFactors = FALSE,
                            data.table = TRUE)

# total[, `:=`(exp_out = paste0(exposure, "_", outcome), effect = "total")]
# total[, mediator := ifelse(exposure != "t2d", "t2d", NA)]
# total[, analysis := paste0(exposure, "_", mediator, "_", outcome)]
# total <- total[exp_out %in% df$exp_out & method %in% c("Wald ratio","Inverse variance weighted")]
# total <- total[, c("exposure","mediator","outcome","b","se","effect","analysis","exp_out")]

total$exp_out <- paste0(total$exposure,"_",total$outcome)
total$effect <- "total"
total$mediator <- ifelse(total$exposure!="t2d","t2d",NA)
total$analysis <- paste0(total$exposure,"_",total$mediator,"_",total$outcome)

total <- total[total$exp_out %in% df$exp_out & total$method %in% c("Wald ratio","Inverse variance weighted"), 
               c("exposure","mediator","outcome","b","se","effect","analysis","exp_out")]

colnames(total) <- colnames(df)

df <- rbind(df,total)
df$exp_out <- NULL

data.table::fwrite(df,"output/twostep_results.csv")