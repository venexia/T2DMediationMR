rm(list=ls())
graphics.off()

df <-  data.table::fread("output/mvmr_results.csv",
                         stringsAsFactors = FALSE,
                         select = c("analysis","exposure","outcome","effect","estimate","se"),
                         data.table = FALSE)

# Multivariable MR with difference of coefficients

df1 <- data.table::fread("output/mvmr_results.csv",
                        stringsAsFactors = FALSE,
                        select = c("analysis","exposure","outcome","effect","estimate","se"),
                        data.table = FALSE)

df1 <- df1[df1$exposure!="t2d",]
df1 <- df1[df1$outcome!="t2d",]

df1 <- tidyr::pivot_wider(df1, names_from = c("effect"), values_from = c("estimate","se"))

df1$estimate <- df1$estimate_total - df1$estimate_direct
df1$se <- sqrt((df1$se_total)^2 + (df1$se_direct)^2)

df1$effect <- "indirect_MVMR"
df1 <- df1[,c("analysis","exposure","outcome","effect","estimate","se")]
colnames(df1) <- colnames(df)

df <- rbind(df,df1)

# Two-step MR with MVMR for second step

df2 <- data.table::fread("output/mvmr_results.csv",
                        stringsAsFactors = FALSE,
                        select = c("analysis","exposure","outcome","effect","estimate","se"),
                        data.table = FALSE)

df2 <- df2[(df2$exposure=="t2d" & df2$effect=="direct") | (df2$outcome=="t2d" & df2$effect=="total"),]
df2$effect <- ifelse(df2$exposure=="t2d","B","A")
df2[,c("exposure","outcome")] <- NULL

df2 <- tidyr::pivot_wider(df2, names_from = c("effect"), values_from = c("estimate","se"))

df2$estimate <- df2$estimate_A * df2$estimate_B
df2$se <- sqrt((df2$se_A)^2 + (df2$se_B)^2)

df2$effect <- "indirect_2stepMVMR"
df2$exposure <- gsub("_.*","",df2$analysis)
df2$outcome <- gsub(".*_","",df2$analysis)
df2<- df2[,c("analysis","exposure","outcome","effect","estimate","se")]
colnames(df2) <- colnames(df)

df <- rbind(df,df2)

# Two-step MR with UVMR for second step

df3 <- data.table::fread("output/mvmr_results.csv",
                         stringsAsFactors = FALSE,
                         select = c("analysis","exposure","outcome","effect","estimate","se"),
                         data.table = FALSE)

df3 <- df3[(df3$exposure=="t2d" & df3$effect=="total") | (df3$outcome=="t2d" & df3$effect=="total"),]
df3 <- df3[df3$outcome %in% c("t2d","pad","cad"),]
df3$effect <- ifelse(df3$exposure=="t2d","B","A")
df3[,c("exposure","outcome")] <- NULL

df3 <- tidyr::pivot_wider(df3, names_from = c("effect"), values_from = c("estimate","se"))

df3$estimate <- df3$estimate_A * df3$estimate_B
df3$se <- sqrt((df3$se_A)^2 + (df3$se_B)^2)

df3$effect <- "indirect_2stepUVMR"
df3$exposure <- gsub("_.*","",df3$analysis)
df3$outcome <- gsub(".*_","",df3$analysis)
df3<- df3[,c("analysis","exposure","outcome","effect","estimate","se")]
colnames(df3) <- colnames(df)

df <- rbind(df,df3)

# Save all estimates

data.table::fwrite(df,"output/mvmr_results_indirect.csv")