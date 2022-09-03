# Gene Signature Application to TCGA
# Copyright (C) ydavidchen 2022

rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("DataLoader.R")
source("GeneSigUtils.R")
set.seed(SEED)

KEY <- "SAMPLE_ID"
load(paste0(DIR_GENES,"enet_object.RData"))

# ------------------------- Subset Expression -------------------------
## Restrict clinical data to ER+/HER2-/LumA IDCs:
cbio <- load_koboldt(censor_mo=MO_CENSOR)
cbio <- subset(cbio, Subtype=="LumA" & 
                     ER_PATH_CALL=="Positive" & 
                     HER2_PATH_CALL=="Negative" & 
                     Oncotree.Code=="IDC")

## Gene expression:
expr <- load_koboldt_expr()
expr$gene_id <- gsub("\\|.*","",expr$gene_id)
expr <- subset(expr, gene_id %in% coefs$Gene)
rownames(expr) <- expr$gene_id
expr$gene_id <- NULL
expr <- as.matrix(expr)
hist(expr)

## Normalize expression data:
expr <- normalize_expr_mat(expr) #right-tailed
hist(expr)

## Apply constraints:
summary(as.numeric(expr))
expr <- custom_winsorize(expr, -WINSOR, WINSOR)
hist(expr)

## Restrict expression data: 
expr <- expr[ , substr(colnames(expr),1,15) %in% cbio[[KEY]] ]

all(cbio[[KEY]] %in% substr(colnames(expr),1,15))
# write.csv(expr, paste0(DIR_OUT, "tcga_gene_sig.csv"), row.names=TRUE, quote=FALSE)

# ------------------------- Cluster Expression -------------------------
res_cl <- samp_annot <- custom_hier_clust(t(expr), CL_PARAMS[1], CL_PARAMS[2])
res_cl[KEY] <- substr(rownames(res_cl), 1, 15)
rownames(res_cl) <- NULL

table(res_cl$Cluster)

wrapper_heatmap(expr, CL_PARAMS[1], CL_PARAMS[2], col_annot=samp_annot)

cbio <- merge(cbio, res_cl, by=KEY)
# write.csv(cbio, paste0(DIR_OUT,"tcga_cluster_info.csv"), row.names=FALSE ,quote=FALSE)
