rm(list=setdiff(ls(), keep))
graphics.off()

# Load features data -----------------------------------------------------------

features <- data.table::fread("data/features_sources.csv",
                              stringsAsFactors = FALSE,
                              data.table = FALSE)

features <- features[features$rsid==TRUE,]

# Load type 2 diabetes GWAS ----------------------------------------------------

t2d <- data.table::fread("data/t2d_rsid.txt",
                         stringsAsFactors = FALSE,
                         data.table = FALSE)

# Load type 2 diabetes instrument ----------------------------------------------

t2d_ins <- data.table::fread("data/instrument-t2d.txt",
                             stringsAsFactors = FALSE,
                             data.table = FALSE)

# Perform analysis for each outcome --------------------------------------------

out = "cad"

# Load outcome GWAS ------------------------------------------------------------

outcome <- data.table::fread(paste0("data/",out,".txt"),
                             stringsAsFactors = FALSE,
                             data.table = FALSE)

# Create results data frame ----------------------------------------------------

results <- data.frame(analysis = character(),
                      direct = logical(),
                      feature = character(),
                      exposure = character(),
                      outcome = character(),
                      beta = numeric(),
                      se = numeric(),
                      pval = numeric(),
                      stringsAsFactors = FALSE)

# Perform analysis for each feature --------------------------------------------

i = 1

## Load feature GWAS -----------------------------------------------------------

feature <- data.table::fread(paste0(path_features_final,features$trait[i],".txt"),
                             stringsAsFactors = FALSE,
                             data.table = FALSE)

## Define feature instrument ---------------------------------------------------

feature_ins <- feature[feature$pval < 5e-8 & !is.na(feature$pval),]
feature_ins$rsid <- feature_ins$SNP
feature_ins <- ieugwasr::ld_clump(dat = feature_ins, clump_kb = 10000, clump_r2 = 0.001)

## Create SNP list -------------------------------------------------------------

snplist <- unique(c(feature_ins$SNP,t2d_ins$SNP))
snplist <- snplist[snplist %in% feature$SNP & snplist %in% t2d$SNP & snplist %in% outcome$SNP]

## Clump SNP list --------------------------------------------------------------

snplist_clumped <- t2d[t2d$SNP %in% snplist,c("SNP","pval")]
colnames(snplist_clumped) <- c("rsid","pval")
snplist_clumped <- ieugwasr::ld_clump(dat = snplist_clumped, clump_kb = 10000, clump_r2 = 0.001)

## Restrict GWAS to relevant SNPs only -----------------------------------------

df <- outcome[outcome$SNP %in% snplist_clumped$rsid,c("SNP","effect_allele","other_allele","beta","se")]
colnames(df) <- c("SNP","effect_allele.outcome","other_allele.outcome","beta.outcome.orig","se.outcome")

df <- merge(df,feature[,c("SNP","effect_allele","other_allele","beta","se")], by = c("SNP"), all.x = TRUE)
colnames(df) <- c("SNP","effect_allele.outcome","other_allele.outcome","beta.outcome.orig","se.outcome",
                  "effect_allele.feature","other_allele.feature","beta.feature.orig","se.feature")

df <- merge(df,t2d[,c("SNP","effect_allele","other_allele","beta","se")], by = c("SNP"), all.x = TRUE)
colnames(df) <- c("SNP","effect_allele.outcome","other_allele.outcome","beta.outcome.orig","se.outcome",
                  "effect_allele.feature","other_allele.feature","beta.feature.orig","se.feature",
                  "effect_allele","other_allele","beta.t2d","se.t2d")

## Harmonise direction of effect -----------------------------------------------

df$beta.feature <- NA
df$beta.feature <- ifelse(df$effect_allele==df$effect_allele.feature & df$other_allele==df$other_allele.feature,df$beta.feature.orig,df$beta.feature)
df$beta.feature <- ifelse(df$effect_allele==df$other_allele.feature & df$other_allele==df$effect_allele.feature,-1*df$beta.feature.orig,df$beta.feature)

df$beta.outcome <- NA
df$beta.outcome <- ifelse(df$effect_allele==df$effect_allele.outcome & df$other_allele==df$other_allele.outcome,df$beta.outcome.orig,df$beta.outcome)
df$beta.outcome <- ifelse(df$effect_allele==df$other_allele.outcome & df$other_allele==df$effect_allele.outcome,-1*df$beta.outcome.orig,df$beta.outcome)

df <- df[,c("SNP","effect_allele","other_allele","beta.t2d","se.t2d","beta.feature","se.feature","beta.outcome","se.outcome")]

## Perform MVMR ----------------------------------------------------------------

### Format data ----------------------------------------------------------------

mvmr_in <- MVMR::format_mvmr(BXGs = df[,c("beta.t2d","beta.feature")],
                            BYG = df[,c("beta.outcome")],
                            seBXGs = df[,c("se.t2d","se.feature")],
                            seBYG = df[,c("se.outcome")],
                            RSID = df[,c("SNP")])

### Test for weak instruments --------------------------------------------------

# strength <- MVMR::strength_mvmr(r_input = mvmr_in, gencov = 0)

### Test for horizontal pleiotropy ---------------------------------------------

# pleiotropy <- MVMR::pleiotropy_mvmr(r_input = mvmr_in, gencov = 0)

### Estimate causal effects ----------------------------------------------------

# mvmr_out <- MVMR::ivw_mvmr(r_inpur = mvmr_in)
mvmr_out <- MVMR::mvmr(r_input = mvmr_in, gencov = 0, weights = 1)

### Record results -------------------------------------------------------------

results[nrow(results)+1,] <- c("MVMR",TRUE,feature$trait[i],"type 2 diabetes",out,mvmr_out$coef[1,c(1,2,4)])
results[nrow(results)+1,] <- c("MVMR",TRUE,feature$trait[i],feature$trait[i],out,mvmr_out$coef[2,c(1,2,4)])

## Perform univariate MR for feature -------------------------------------------

### Format data ----------------------------------------------------------------

uvmr_in <- MendelianRandomization::mr_input(bx = df[df$SNP %in% feature_ins$SNP,]$beta.feature,
                                           bxse = df[df$SNP %in% feature_ins$SNP,]$se.feature,
                                           by = df[df$SNP %in% feature_ins$SNP,]$beta.outcome,
                                           byse = df[df$SNP %in% feature_ins$SNP,]$se.outcome)

### Estimate causal effects ----------------------------------------------------

uvmr_out <- MendelianRandomization::mr_ivw(uvmr_in)

### Record results -------------------------------------------------------------

results[nrow(results)+1,] <- c("UVMR",NA,feature$trait[i],feature$trait[i],out,uvmr_out@Estimate,uvmr_out@StdError,uvmr_out@Pvalue)

## Perform univariate MR for type 2 diabetes -----------------------------------

### Format data ----------------------------------------------------------------

uvmr_in <- MendelianRandomization::mr_input(bx = df[df$SNP %in% t2d_ins$SNP,]$beta.t2d,
                                            bxse = df[df$SNP %in% t2d_ins$SNP,]$se.t2d,
                                            by = df[df$SNP %in% t2d_ins$SNP,]$beta.outcome,
                                            byse = df[df$SNP %in% t2d_ins$SNP,]$se.outcome)

### Estimate causal effects ----------------------------------------------------

uvmr_out <- MendelianRandomization::mr_ivw(uvmr_in)

### Record results -------------------------------------------------------------

results[nrow(results)+1,] <- c("UVMR",NA,feature$trait[i],"type 2 diabetes",out,uvmr_out@Estimate,uvmr_out@StdError,uvmr_out@Pvalue)

## Calculate indirect effects --------------------------------------------------

direct_beta.t2d <- results[results$analysis=="MVMR" & results$feature==feature$trait[i] & results$exposure=="type 2 diabetes",]$beta
univar_beta.t2d <- results[results$analysis=="UVMR" & results$feature==feature$trait[i] & results$exposure=="type 2 diabetes",]$beta

direct_se.t2d <- results[results$analysis=="MVMR" & results$feature==feature$trait[i] & results$exposure=="type 2 diabetes",]$se
univar_se.t2d <- results[results$analysis=="UVMR" & results$feature==feature$trait[i] & results$exposure=="type 2 diabetes",]$se

indirect_beta.t2d <- univar_beta.t2d - direct_beta.t2d 
indirect_se.t2d <- sqrt((uv[uv$exposure==i,]$se)^2 + (res[res$exposure==i,]$se)^2)
results[nrow(results)+1,] <- c("MVMR",FALSE,feature$trait[i],"type 2 diabetes",out,indirect_beta.t2d,indirect_se.t2d,NA)

indirect_beta.feature <-
indirect_se.feature <- 
results[nrow(results)+1,] <- c("MVMR",FALSE,feature$trait[i],feature$trait[i],out,indirect_beta.feature,indirect_se.feature,NA)
