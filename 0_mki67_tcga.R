# MKI67 in TCGA LumA vs. LumB
# Copyright (C) ydavidchen 2022

rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("DataLoader.R")
source("GeneSigUtils.R")
LUMS <- c("LumA", "LumB")
S_KEY <- "SAMPLE_ID"

## Normalized MKI67 RNA-seq expression from cBioPortal:
mki67 <- load_cbio_query(paste0(TCGA_DIR,"query1_TCGABRCA/mRNA expression z-scores relative to all samples (log RNA Seq V2 RSEM).txt"), "mRNA")

## Restrict to ER+/HER2-/Luminal IDCs:
cbio <- load_koboldt()
cbio <- subset(cbio, Subtype %in% LUMS & 
                 ER_PATH_CALL=="Positive" & 
                 HER2_PATH_CALL=="Negative" & 
                 Oncotree.Code=="IDC")

cbio <- merge(cbio, mki67, by=S_KEY)

## Comparisons:
ggplot(cbio, aes(Subtype, mRNA)) + 
  geom_boxplot(outlier.shape=NA) + 
  geom_jitter(width=0.2) +
  labs(y="Normalized MKI67", title="TCGA") +
  BOXPLOT_THEME

t.test(
  cbio$mRNA[cbio$Subtype==LUMS[1]],
  cbio$mRNA[cbio$Subtype==LUMS[2]]
)
