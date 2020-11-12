tmp <- df[df$method %in% c("Inverse variance weighted","Wald ratio"), c("exposure","outcome","evidence")]

tmp <- tmp[tmp$exposure %in% features | tmp$outcome %in% features,]

tmp$feature <- NA
tmp$feature <- ifelse(tmp$outcome %in% c("t2d_udler","t2d","cad","pad"),tmp$exposure,tmp$feature)
tmp$feature <- ifelse(tmp$exposure %in% c("t2d_udler","t2d","cad","pad"),tmp$outcome,tmp$feature)

tmp$exposure <- ifelse(tmp$exposure==tmp$feature, "feature", tmp$exposure)
tmp$outcome <- ifelse(tmp$outcome==tmp$feature, "feature", tmp$outcome)
tmp$analysis <- paste(tmp$exposure, tmp$outcome, sep = "_")

tmp <- tmp[,c("analysis","feature","evidence")]

tmp <- tidyr::pivot_wider(tmp, names_from = "analysis", values_from = "evidence")
tmp$bidir <- tmp$feature_t2d==TRUE & tmp$t2d_feature==TRUE