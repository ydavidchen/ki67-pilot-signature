# Gene Signature: Discovery Analysis
# Copyright (C) ydavidchen 2022

rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("GeneSigUtils.R")
library(glmnet)
set.seed(SEED)

## Constants & data:
LOW_GROUP <- "********** MASKED **********"
HIGH_GROUP <- "********** MASKED **********"

wu <- read.csv(paste0(DIR_WU,"GSE176078_Wu_etal_2021_bulkRNAseq_raw_counts_resaved.csv"))
commonSet <- read.table(PATH_COMMONSET)$V1

## 1. Subset to gene universe shared across cohorts:
wu <- subset(wu, wus95idmgmtGenes %in% commonSet)
rownames(wu) <- wu$wus95idmgmtGenes
wu$wus95idmgmtGenes <- NULL

## 2a. Sample-wise normalization (across ALL samples)
hist( as.numeric(as.matrix(wu)) ) #counts

wu <- normalize_expr_mat(as.matrix(wu))
hist(wu)

## 2b. Winsorize normalized data:
summary(as.numeric(wu))
wu <- custom_winsorize(wu, -WINSOR, WINSOR)
hist(wu)

## 3. Code regression target;
wu <- t(wu)
y <- factor( c(rep(1, length(HIGH_GROUP)), rep(0, length(LOW_GROUP))) )
wu <- wu[c(HIGH_GROUP,LOW_GROUP), ] #reorder

## 4. Build model
stopifnot(identical(rownames(wu), c(HIGH_GROUP,LOW_GROUP))) #checkpoint

eNet <- glmnet(wu, y, family="binomial", alpha=0.5)
plot(eNet)
lambda_best <- min(eNet$lambda)
lambda_best #0.009797275

## 5. Extract coefficients:
coefs <- coef(eNet, s=lambda_best)
coefs <- data.frame(
  Gene = coefs@Dimnames[[1]][1+coefs@i],
  Coef = coefs@x
)

## Visualization
samp_annot <- data.frame(
  row.names = c(HIGH_GROUP, LOW_GROUP),
  HighKi67 = ifelse(y==1, "Yes", "No")
)

wrapper_heatmap(
  t(wu[,colnames(wu) %in% coefs$Gene]), 
  clustering_distance_rows = CL_PARAMS[1],
  clustering_method = CL_PARAMS[2],
  col_annot = samp_annot, 
  showcn = TRUE
)

ggplot(aes(reorder(Gene, Coef), Coef), data=coefs) +
  geom_bar(stat="identity") +
  labs(x="Gene", y="Elastic-net coefficient") +
  scale_y_continuous(breaks=seq(-1.3,0.5,0.2)) +
  BARPLOT_THEME

## Export matrix & object:
# save(list=c("coefs","commonSet","eNet"), file=paste0(DIR_OUT,"enet_object.RData"), compress=TRUE)
# write.csv(coefs, paste0(DIR_OUT,"enet_coefs.csv"), row.names=FALSE, quote=FALSE)
