# Gene Signature Application to METABRIC
# Copyright (C) ydavidchen 2022

rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("DataLoader.R")
source("GeneSigUtils.R")

KEY <- "PATIENT_ID"
load(paste0(DIR_GENES,"enet_object.RData"))

# ------------------------- Subset Expression -------------------------
## Restrict clinical data to ER+/HER2-/LumA IDCs:
cbio <- load_curtis(censor_months=MO_CENSOR)
cbio <- subset(cbio, CLAUDIN_SUBTYPE=="LumA" & 
                     ER_IHC=="Positve" & 
                     HER2_STATUS=="Negative" & 
                     HISTOLOGICAL_SUBTYPE=="Ductal/NST") 

## Load expression & scale sample-wise manually:
expr <- load_curtis_expr(zScoreFile=FALSE) #transcript count
expr <- subset(expr, Hugo_Symbol %in% coefs$Gene)
sum(duplicated(expr$Hugo_Symbol)) #2
expr <- aggregate(. ~  Hugo_Symbol, data=expr, FUN=mean)
rownames(expr) <- expr$Hugo_Symbol
expr$Hugo_Symbol <- NULL
expr <- as.matrix(expr)
hist(expr)

## Sample-wise normalize:
expr <- normalize_expr_mat(expr)
hist(expr)

## Apply constraints:
expr <- custom_winsorize(expr, -WINSOR, WINSOR)
hist(expr)

## Restrict expression data / mutually subset:
expr <- expr[ , colnames(expr) %in% cbio[[KEY]] ]
cbio <- subset(cbio, cbio[[KEY]] %in% colnames(expr))
# write.csv(expr, paste0(DIR_OUT, "metabric_gene_sig.csv"), row.names=TRUE, quote=FALSE)

# ------------------------- Cluster Expression -------------------------
res_cl <- samp_annot <- custom_hier_clust(t(expr), CL_PARAMS[1], CL_PARAMS[2])
res_cl[KEY] <- gsub(".", "-", rownames(res_cl), fixed=TRUE)
rownames(res_cl) <- NULL

table(res_cl$Cluster)

wrapper_heatmap(expr, CL_PARAMS[1], CL_PARAMS[2], col_annot=samp_annot)

cbio <- merge(cbio, res_cl, by=KEY)
# write.csv(cbio, paste0(DIR_OUT,"metabric_cluster_info.csv"), row.names=FALSE ,quote=FALSE)
