uvmr <- data.table::fread("output/uvmr_results.csv", data.table = FALSE)
uvmr$analysis <- paste(uvmr$exposure,uvmr$outcome,sep="_")
uvmr <- uvmr[uvmr$method %in% c("Inverse variance weighted","Wald ratio") & 
               uvmr$outcome!="t2d" &
               !(uvmr$exposure %in% c("t2d","t2d_udler")),]
uvmr$effect <- "total"
uvmr$estimate <- uvmr$b

mvmr <- data.table::fread("output/mvmr_results_indirect.csv", data.table = FALSE)
mvmr <- mvmr[mvmr$outcome!="t2d" & !(mvmr$exposure %in% c("t2d","t2d_udler")),]

df <- rbind(uvmr[,c("analysis","exposure","outcome","effect","estimate","se")],
            mvmr[,c("analysis","exposure","outcome","effect","estimate","se")])


df <- tidyr::pivot_wider(df, names_from = c("effect"), values_from = c("estimate","se"))

df$analysis <- NULL
df <- df[,c("exposure","outcome",
            "estimate_total","se_total",
            "estimate_total_restricted","se_total_restricted",
            "estimate_direct","se_direct",
            "estimate_indirect_MVMR","se_indirect_MVMR",
            "estimate_indirect_2stepMVMR","se_indirect_2stepMVMR",
            "estimate_indirect_2stepUVMR","se_indirect_2stepUVMR")]