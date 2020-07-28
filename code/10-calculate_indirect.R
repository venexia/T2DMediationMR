rm(list=ls())
graphics.off()

df <-  data.table::fread("output/mvmr_results.csv",
                         stringsAsFactors = FALSE,
                         data.table = FALSE)

# Multivariable MR with difference of coefficients

df1 <- df[,c("analysis","exposure","outcome","effect","estimate","se")]
df1 <- df1[df1$exposure!="t2d",]
df1 <- df1[df1$outcome!="t2d",]

df1 <- tidyr::pivot_wider(df1, names_from = c("effect"), values_from = c("estimate","se"))

df1$estimate <- df1$estimate_total_restricted - df1$estimate_direct
df1$se <- sqrt((df1$se_total_restricted)^2 + (df1$se_direct)^2)

df1$effect <- "indirect_MVMR"
df1 <- df1[,c("analysis","exposure","outcome","effect","estimate","se")]
colnames(df1) <- colnames(df)

df <- plyr::rbind.fill(df,df1)

# Two-step MR with MVMR for second step

for (i in c("cad","pad")) {
  
  df2 <- df[grepl(i,df$analysis),c("analysis","exposure","outcome","effect","estimate","se")]
  df2 <- df2[!(df2$exposure %in% c("cad","pad")),]
  df2 <- df2[(df2$exposure=="t2d" & df2$effect=="direct") | (df2$outcome=="t2d" & df2$effect=="total_restricted"),]
  df2 <- df2[df2$outcome %in% c("t2d",i),]
  df2$effect <- ifelse(df2$exposure=="t2d","B","A")
  df2$analysis <- paste(substr(df2$analysis,1,nchar(df2$analysis)-4), i, sep = "_")
  df2[,c("exposure","outcome")] <- NULL
  
  df2 <- tidyr::pivot_wider(df2, names_from = c("effect"), values_from = c("estimate","se"))
  
  df2$estimate <- df2$estimate_A * df2$estimate_B
  df2$se <- sqrt((df2$se_A)^2 + (df2$se_B)^2)
  
  df2$effect <- "indirect_2stepMVMR"
  df2$exposure <- gsub("_.*","",df2$analysis)
  df2$outcome <- gsub(".*_","",df2$analysis)
  df2<- df2[,c("analysis","exposure","outcome","effect","estimate","se")]
  colnames(df2) <- colnames(df)
  
  df <- plyr::rbind.fill(df,df2)
  
}

# Two-step MR with UVMR for second step


for (i in c("cad","pad")) {

  df3 <- data.table::fread("output/uvmr_results.csv",
                           stringsAsFactors = FALSE,
                           select = c("method","exposure","outcome","b","se"),
                           data.table = FALSE)
  
  colnames(df3) <- c("method","exposure","outcome","estimate.uvmr","se.uvmr")
  
  estimate.B <- df3[df3$exposure=="t2d" & df3$outcome==i & df3$method %in% c("Inverse variance weighted","Wald ratio"),]$estimate.uvmr
  se.B <- df3[df3$exposure=="t2d" & df3$outcome==i & df3$method %in% c("Inverse variance weighted","Wald ratio"),]$se.uvmr
  
  df3 <- df3[df3$method %in% c("Inverse variance weighted","Wald ratio") &
               df3$outcome==i &
               !(df3$exposure %in% c("t2d","t2d_udler")),]
  
  df3$estimate <- df3$estimate.uvmr * estimate.B
  df3$se <- sqrt((df3$se.uvmr)^2 + (se.B)^2)
  
  df3$effect <- "indirect_2stepUVMR"
  df3$analysis <- paste(df3$exposure, df3$outcome, sep="_")
  
  df3 <- df3[,c("analysis","exposure","outcome","effect","estimate","se")]
  
  df <- plyr::rbind.fill(df,df3)

}

# Save all estimates

data.table::fwrite(df,"output/mvmr_results_indirect.csv")