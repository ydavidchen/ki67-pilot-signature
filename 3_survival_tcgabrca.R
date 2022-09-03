# Clinical & Survival in TCGA
# Copyright (C) ydavidchen 2022

rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("GeneSigUtils.R")
library(tableone)
library(survival)
library(survminer)

df_surv <- read.csv(paste0(DIR_GENES,"tcga_cluster_info.csv"))

#------------------------------ Clinical ------------------------------
nVars <- c("Diagnosis.Age", "TMB..nonsynonymous.")
fVars <- c("menopause_status", "BinaryAJCC", "Race.Category")

t1 <- CreateTableOne(
  vars = c(nVars, fVars),
  strata = "Cluster",
  data = df_surv,
  factorVars = fVars,
  includeNA = TRUE,
  test = TRUE
)
t1 <- print(t1, showAllLevels=TRUE)
# write.csv(t1, paste0(DIR_OUT,"tcga_tableone.csv"))

#------------------------------ Kaplan-Meier ------------------------------
osFit <- survfit(Surv(OS_time, OS_status)~Cluster, data=df_surv)
ggsurvplot(osFit, df_surv, palette=SURV_COLS, ggtheme=SURV_THEME, pval=TRUE)

dfsFit <- survfit(Surv(DFS_time, DFS_status)~Cluster, data=df_surv)
ggsurvplot(dfsFit, df_surv, palette=SURV_COLS, ggtheme=SURV_THEME, pval=TRUE)

#------------------------------ CoxPH ------------------------------
df_surv$Cluster <- factor(df_surv$Cluster, levels=c("Cluster2","Cluster1"))

osCox <- coxph(Surv(OS_time, OS_status)~Cluster+Diagnosis.Age+BinaryAJCC, data=df_surv)
summary(osCox)

dfsCox <- coxph(Surv(DFS_time, DFS_status)~Cluster+Diagnosis.Age+BinaryAJCC, data=df_surv)
summary(dfsCox)
