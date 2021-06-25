rm(list=ls())
graphics.off()

# Load mvmr analysis summary ---------------------------------------------------

df <- data.table::fread("data/mvmr_summary.csv", data.table = FALSE)

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
df <- merge(df, mvmr, by = c("exposure","mediator","outcome"), all.x = TRUE)

# Prepare UVMR results ---------------------------------------------------------

uvmr_fwd <-  data.table::fread("output/results_fwd.csv",
                                select = c("exposure","outcome","method","b","se"),
                                data.table = FALSE)

uvmr_bkwd <-  data.table::fread("output/results_bkwd.csv",
                           select = c("exposure","outcome","method","b","se"),
                           data.table = FALSE)

uvmr <- rbind(uvmr_fwd, uvmr_bkwd)

uvmr <- uvmr[uvmr$method %in% c("Inverse variance weighted","Wald ratio"), c("exposure","outcome","b","se")]

colnames(uvmr) <- c("exposure","mediator","estimate_A","se_A")

df <- merge(df, uvmr, by = c("exposure","mediator"), all.x = TRUE)

# Restrict to non-missing estimates --------------------------------------------

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
df$Qstat <- NA
df$Qpval <- NA
df$Qdf <- NA
df$condF <- NA
  
# Add direct effect results ----------------------------------------------------

direct <-  data.table::fread("output/mvmr_results.csv",
                             select = c("analysis","exposure","outcome","effect","estimate","se","Qstat","Qpval","Qdf","condF"),
                             stringsAsFactors = FALSE,
                             data.table = FALSE)

direct <- direct[direct$effect=="direct",]
direct <- tidyr::separate(direct, analysis, sep = "/", into = c("exposure1","exposure2","outcome1"))
direct$mediator <- ""
direct$mediator <- ifelse(direct$exposure1==direct$exposure,direct$exposure2,direct$mediator)
direct$mediator <- ifelse(direct$exposure2==direct$exposure,direct$exposure1,direct$mediator)

direct <- direct[,c("exposure","mediator","outcome","estimate","se","effect","Qstat","Qpval","Qdf","condF")]

direct <- merge(df[,c("exposure","mediator","outcome")], direct, by = c("exposure","mediator","outcome"), all.x = TRUE)
df <- rbind(df,direct)

# Add total effect results -----------------------------------------------------

total_fwd <-  data.table::fread("output/results_fwd.csv",
                               select = c("exposure","outcome","method","b","se"),
                               data.table = FALSE)

total_bkwd <-  data.table::fread("output/results_bkwd.csv",
                                select = c("exposure","outcome","method","b","se"),
                                data.table = FALSE)

total <- rbind(total_fwd, total_bkwd)

total <- total[total$method %in% c("Wald ratio","Inverse variance weighted"), c("exposure","outcome","b","se")]

colnames(total) <- c("exposure","outcome","estimate","se")

total$effect <- "total"

total <- merge(unique(df[,c("exposure","mediator","outcome")]), total, 
               by = c("exposure","outcome"), all.x = TRUE)

total$Qstat <- NA
total$Qpval <- NA
total$Qdf <- NA
total$condF <- NA

total <- total[,colnames(df)]

df <- rbind(df,total)

data.table::fwrite(df,"output/twostep_results.csv")