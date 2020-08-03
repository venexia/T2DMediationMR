rm(list=ls())
graphics.off()

# Prepare MVMR results ---------------------------------------------------------

mvmr <-  data.table::fread("output/mvmr_results.csv",
                           select = c("analysis","exposure","outcome","effect","estimate","se"),
                           stringsAsFactors = FALSE,
                           data.table = FALSE)

mvmr <- mvmr[mvmr$outcome %in% c("pad","cad") & mvmr$effect=="direct" & mvmr$exposure=="t2d",]
mvmr$effect <- NULL
colnames(mvmr) <- c("analysis","mediator","outcome","estimate_B","se_B")

mvmr$exposure <- substr(mvmr$analysis,1,nchar(mvmr$analysis)-4)
mvmr$analysis <- NULL

# Prepare UVMR results ---------------------------------------------------------

uvmr <-  data.table::fread("output/uvmr_results.csv",
                           select = c("exposure","outcome","method","b","se"),
                           stringsAsFactors = FALSE,
                           data.table = FALSE)

uvmr <- uvmr[uvmr$method %in% c("Inverse variance weighted","Wald ratio") & uvmr$outcome=="t2d",c("exposure","outcome","b","se")]

colnames(uvmr) <- c("exposure","mediator","estimate_A","se_A")

uvmr <- uvmr[uvmr$exposure %in% mvmr$exposure,]

# Combine estimates ------------------------------------------------------------

df <- merge(uvmr, mvmr, by = c("exposure","mediator"), all = TRUE)

# Calculate indirect effects ---------------------------------------------------

df$estimate <- df$estimate_A * df$estimate_B
# df$se <- abs(df$estimate) * sqrt((df$se_A/df$estimate_A)^2 + (df$se_B/df$estimate_B)^2)
df$se <- sqrt((df$se_A)^2 + (df$se_B)^2)

df <- df[,c("exposure","mediator","outcome","estimate_A","se_A","estimate_B","se_B","estimate","se")]
df$pc_A <- df$se_A/abs(df$estimate_A)
df$pc_B <- df$se_B/abs(df$estimate_B)
df$pc <- df$se/abs(df$estimate)

df <- df[,c("exposure","outcome","estimate","se")]
df$analysis <- paste0(df$exposure,"_",df$outcome)
df$effect <- "indirect"

data.table::fwrite(df,"output/twostep_results.csv")

# Save all estimates -----------------------------------------------------------

mvmr <-  data.table::fread("output/mvmr_results.csv",
                           stringsAsFactors = FALSE,
                           data.table = FALSE)

uvmr <-  data.table::fread("output/uvmr_results.csv",
                           stringsAsFactors = FALSE,
                           data.table = FALSE)

uvmr <- uvmr[uvmr$method %in% c("Inverse variance weighted","Wald ratio"),]
uvmr$analysis <- paste0(uvmr$exposure,"_",uvmr$outcome)
uvmr$effect <- "total"
uvmr$estimate <- uvmr$b
uvmr$snps <- uvmr$snps_in
uvmr[,c("id.exposure","id.outcome","b","snps_in","snps_mr","method")] <- NULL

all <- plyr::rbind.fill(uvmr,mvmr)

all <- plyr::rbind.fill(all,df)

data.table::fwrite(all,"output/all_results.csv")