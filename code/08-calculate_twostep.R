rm(list=ls())
graphics.off()

# Load filtered features list --------------------------------------------------

features <- data.table::fread("output/evidence_summary.csv", data.table = FALSE)

features$exposure <- NA
features$exposure <- ifelse(features$feature_t2d_evidence==TRUE & features$t2d_feature_evidence==FALSE,features$trait,features$exposure)
features$exposure <- ifelse(features$feature_t2d_evidence==FALSE & features$t2d_feature_evidence==TRUE,"t2d",features$exposure)

features$mediator <- NA
features$mediator <- ifelse(features$feature_t2d_evidence==TRUE & features$t2d_feature_evidence==FALSE,"t2d",features$mediator)
features$mediator <- ifelse(features$feature_t2d_evidence==FALSE & features$t2d_feature_evidence==TRUE,features$trait,features$mediator)

features <- tidyr::pivot_longer(features,
                                cols = c("feature_pad_evidence","feature_cad_evidence"),
                                names_to = "analysis",
                                values_to = "evidence")

features$outcome <- NA
features$outcome <- ifelse(features$analysis=="feature_pad_evidence" & features$evidence==TRUE,"pad",features$outcome)
features$outcome <- ifelse(features$analysis=="feature_cad_evidence" & features$evidence==TRUE,"cad",features$outcome)

features <- features[,c("exposure","mediator","outcome")]
features <- na.omit(features)

tmp <- features[features$mediator=="t2d",]
tmp$mediator <- "t2d_linear"
features <- rbind(features,tmp)

tmp <- features[features$exposure=="t2d",]
tmp$exposure <- "t2d_linear"
features <- rbind(features,tmp)

# Prepare MVMR results ---------------------------------------------------------

mvmr <-  data.table::fread("output/mvmr_results.csv",
                           select = c("analysis","exposure","outcome","effect","estimate","se"),
                           stringsAsFactors = FALSE,
                           data.table = FALSE)

mvmr <- mvmr[mvmr$effect=="direct",]
mvmr$effect <- NULL
colnames(mvmr) <- c("analysis","mediator","outcome","estimate_B","se_B")

mvmr <- tidyr::separate(mvmr, analysis, sep = "/", into = c("exposure1","exposure2","outcome1"))
mvmr$exposure <- ""
mvmr$exposure <- ifelse(mvmr$exposure1==mvmr$mediator,mvmr$exposure2,mvmr$exposure)
mvmr$exposure <- ifelse(mvmr$exposure2==mvmr$mediator,mvmr$exposure1,mvmr$exposure)
mvmr <- mvmr[,c("exposure","mediator","outcome","estimate_B","se_B")]

mvmr <- merge(features,mvmr,all.x = TRUE, by = c("exposure","mediator","outcome"))
mvmr$exposure_mediator <- paste0(mvmr$exposure,"/",mvmr$mediator)

# Prepare UVMR results ---------------------------------------------------------

uvmr <-  data.table::fread("output/results.csv",
                           select = c("exposure","outcome","method","b","se"),
                           stringsAsFactors = FALSE,
                           data.table = FALSE)

uvmr <- uvmr[uvmr$method %in% c("Inverse variance weighted","Wald ratio"), c("exposure","outcome","b","se")]

colnames(uvmr) <- c("exposure","mediator","estimate_A","se_A")
uvmr$exposure_mediator <- paste0(uvmr$exposure,"/",uvmr$mediator)

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
df$analysis <- paste0(df$exposure,"/",df$mediator,"/",df$outcome)
df$Qstat <- NA
df$Qpval <- NA
df$Qdf <- NA
df$condF <- NA
  
# Add direct effect results ----------------------------------------------------

direct <-  data.table::fread("output/mvmr_results.csv",
                             select = c("analysis","exposure","outcome","effect","estimate","se","Qstat","Qpval","Qdf","condF"),
                             stringsAsFactors = FALSE,
                             data.table = FALSE)

direct <- tidyr::separate(direct, analysis, sep = "/", into = c("exposure1","exposure2","outcome1"))
direct$mediator <- ""
direct$mediator <- ifelse(direct$exposure1==direct$exposure,direct$exposure2,direct$mediator)
direct$mediator <- ifelse(direct$exposure2==direct$exposure,direct$exposure1,direct$mediator)
direct$analysis <- paste0(direct$exposure,"/",direct$mediator,"/",direct$outcome)

direct <- direct[direct$analysis %in% df$analysis & direct$effect=="direct", 
                 c("exposure","mediator","outcome","estimate","se","effect","analysis","Qstat","Qpval","Qdf","condF")]

df <- rbind(df,direct)

# Add total effect results -----------------------------------------------------

df$exp_out <- paste0(df$exposure,"/",df$outcome)

total <-  data.table::fread("output/results.csv",
                            stringsAsFactors = FALSE,
                            data.table = TRUE)

total$exp_out <- paste0(total$exposure,"/",total$outcome)
total$effect <- "total"
total$mediator <- ifelse(total$exposure!="t2d","t2d",NA)

total <- total[total$exp_out %in% df$exp_out & total$method %in% c("Wald ratio","Inverse variance weighted"), 
               c("exposure","mediator","outcome","b","se","effect","exp_out")]

tmp <- total
tmp$mediator <- ifelse(tmp$mediator=="t2d","t2d_linear",tmp$mediator)
total <- rbind(total,tmp)

total$analysis <- paste0(total$exposure,"/",total$mediator,"/",total$outcome)
colnames(total) <- c("exposure","mediator","outcome","estimate","se","effect","analysis","exp_out")
total$Qstat <- NA
total$Qpval <- NA
total$Qdf <- NA
total$condF <- NA

df <- rbind(df,total)
df$exp_out <- NULL

data.table::fwrite(df,"output/twostep_results.csv")