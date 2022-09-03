# Clinical & Survival in METABRIC
# Copyright (C) ydavidchen 2022

rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("GeneSigUtils.R")
library(tableone)
library(survival)
library(survminer)

df_surv <- read.csv(paste0(DIR_GENES,"metabric_cluster_info.csv"))

# ------------------------- Clinical -------------------------
nVars <- c("AGE_AT_DIAGNOSIS", "TMB_NONSYNONYMOUS")
fVars <- c("INFERRED_MENOPAUSAL_STATE", "BinaryStage")

t1 <- CreateTableOne(
  vars = c(nVars, fVars),
  strata = "Cluster",
  data = df_surv,
  factorVars = fVars,
  includeNA = TRUE,
  test = TRUE
)
t1 <- print(t1, showAllLevels=TRUE)
t.test(AGE_AT_DIAGNOSIS ~ Cluster, data=df_surv) #exact pval
# write.csv(t1, paste0(DIR_OUT,"./metabric_tableone.csv"))

# ------------------------- Kaplan-Meier -------------------------
osFit <- survfit(Surv(OS_time, OS_status)~Cluster, data=df_surv)
ggsurvplot(osFit, df_surv, palette=SURV_COLS, ggtheme=SURV_THEME, pval=TRUE)

rfsFit <- survfit(Surv(RFS_time, RFS_status)~Cluster, data=df_surv)
ggsurvplot(rfsFit, df_surv, palette=SURV_COLS, ggtheme=SURV_THEME, pval=TRUE)

# ------------------------- CoxPH -------------------------
df_surv$Cluster <- factor(df_surv$Cluster, levels=c("Cluster2","Cluster1"))

osCox <- coxph(Surv(OS_time, OS_status)~Cluster+AGE_AT_DIAGNOSIS+BinaryStage, data=df_surv)
summary(osCox)

rfsCox <- coxph(Surv(RFS_time, RFS_status)~Cluster+AGE_AT_DIAGNOSIS+BinaryStage, data=df_surv)
summary(rfsCox)
