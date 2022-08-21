# MKI67 in METABRIC LumA vs. LumB
# Copyright (C) ydavidchen 2022

rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("DataLoader.R")
source("GeneSigUtils.R")

KEY <- "PATIENT_ID"
LUMS <- c("LumA", "LumB")

## Load MKI67:
mki67 <- load_curtis_expr(zScoreFile=TRUE)
mki67 <- subset(mki67, Hugo_Symbol=="MKI67")
rownames(mki67) <- mki67$Hugo_Symbol
mki67$Entrez_Gene_Id <- mki67$Hugo_Symbol <- NULL
mki67 <- as.data.frame(t(mki67))
mki67[KEY] <- rownames(mki67)

## Restrict to ER+/HER2-/Luminal IDCs:
cbio <- load_curtis()
cbio <- subset(cbio, CLAUDIN_SUBTYPE %in% LUMS & 
                     ER_IHC=="Positve" & 
                     HER2_STATUS=="Negative" & 
                     HISTOLOGICAL_SUBTYPE=="Ductal/NST") 

## Compare MKI67:
cbio <- merge(cbio, mki67, by=KEY)

ggplot(cbio, aes(CLAUDIN_SUBTYPE, MKI67)) +
  geom_boxplot(outlier.shape=NA) + 
  geom_jitter(width=0.2) +
  labs(y="Normalized MKI67", title="METABRIC") +
  BOXPLOT_THEME

t.test(
  cbio$MKI67[cbio$CLAUDIN_SUBTYPE==LUMS[1]],
  cbio$MKI67[cbio$CLAUDIN_SUBTYPE==LUMS[2]]
)
