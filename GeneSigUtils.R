# Utility Module for Gene Signature Analysis
# Copyright (C) ydavidchen 2022

library(matrixStats)
library(ggplot2)
library(pheatmap)

## Path variables (MASKED):
DIR_WU <- "*****MASKED*****"
PATH_COMMONSET <- "*****MASKED*****"
DIR_GENES <- "*****MASKED*****"
DIR_OUT <- "*****MASKED*****"

## Constants for analysis:
WINSOR <- 3.5
CL_PARAMS <- c("euclidean", "ward.D")
SEED <- 42
MO_CENSOR <- 10 * 12

## Helper functions:
normalize_expr_mat <- function(mat) {
  #'@description Sample-wise mean-centering, variance-normalization
  #'@details Reminder: Row=Gene, column=sample
  mu <- as.matrix(rowMeans(mat))
  sd <- as.matrix(rowSds(mat))
  mat <- sweep(mat, 1, mu, FUN="-")
  mat <- sweep(mat, 1, sd, FUN="/")
  return(mat)
}

custom_winsorize <- function(mat, lower=NULL, upper=NULL) {
  #'@description Cap matrix value at lower & upper bounds
  if(!is.null(upper)) {
    mat[mat > upper] <- upper
  }
  if(!is.null(lower)) {
    mat[mat < lower] <- lower
  }
  return(mat)
}

custom_hier_clust <- function(tMat, method_dist, method_hc, num_cl=2) {
  #'@param tMat Expression matrix where row=Samples(observations), col=genes(features)
  myDist <- dist(tMat, method=method_dist)
  myHcl <- hclust(myDist, method=method_hc)
  
  plot(myHcl, labels=FALSE)
  rect.hclust(myHcl, k=num_cl)
  
  membership <- data.frame(cutree(myHcl, k=num_cl))
  colnames(membership) <- "Cluster"
  membership$Cluster <- paste0("Cluster", membership$Cluster)
  return(membership)
}

add_num_to_lab <- function(col, newLine=TRUE) {
  #'@description Creates new column/vector adding n
  tabl <- as.data.frame(table(col)) #S3 method
  tabl$n <- paste0("(n=", tabl$Freq, ")")
  s <- ifelse(newLine, "\n", " ") #space or newline
  tabl$newLab <- paste0(tabl$col, s, tabl$n)
  map <- setNames(tabl$newLab, tabl$col)
  return(map[col])
}

wrapper_heatmap <- function(mat, clustering_distance_rows, clustering_method,
                            zscore=FALSE, row_annot=NA, col_annot=NA, 
                            showrn=FALSE, showcn=FALSE) {
  pheatmap(
    mat,
    show_rownames = showrn,
    show_colnames = showcn,
    clustering_distance_rows=clustering_distance_rows,
    clustering_method=clustering_method,
    border_color = NA,
    fontsize = 12,
    color = HEAT_COLS,
    annotation_colors = ANN_COLORS,
    scale = ifelse(zscore, "row", "none"),
    annotation_col = col_annot,
    annotation_row = row_annot
  )
}

## Constants for data visualization:
HEAT_COLS <- colorRampPalette(c("blue","lightgray","red"))(1024)
SURV_COLS <- c("black","gray")
ANN_COLORS <- list(
  Cluster = c(Cluster1="black", Cluster2="lightgray"),
  HighKi67 = c(Yes="black", No="lightgray")
)

SURV_THEME <- theme_classic() +
  theme(axis.text=element_text(size=20,color="black"), 
        axis.title=element_text(size=20,color="black"),
        title=element_text(size=20,color="black"),
        legend.title=element_text(size=16,color="black",face="bold"), legend.text=element_text(size=16,color="black",face="bold"), 
        legend.box.background=element_rect(colour="black",size=2),
        text=element_text(size=20),
        strip.text.x=element_text(size=20,colour="black",face="bold"))

SCATTER_THEME <- theme_bw() + 
  theme(axis.text.x=element_text(size=20,color="black"), axis.title.x=element_text(size=20,color="black"),
        axis.text.y=element_text(size=20,color="black"), axis.title.y=element_text(size=20,color="black"),
        title = element_text(size=20,color="black"), 
        panel.spacing = unit(0.5, "lines"), panel.border = element_blank(), axis.line=element_line(color="black"),
        strip.text.x=element_text(size=20,color="black",face="bold"), strip.background=element_rect(fill="gray95"),
        legend.position="top", legend.title=element_blank(), legend.text=element_text(size=15,color="black")) 

BOXPLOT_THEME <- theme_bw() +
  theme(axis.text.x=element_text(size=15,color="black"), axis.text.y=element_text(size=15,color="black"),
        axis.title.x=element_blank(), axis.title.y=element_text(size=20,color="black"),
        strip.text.x=element_text(size=20,color="black",face="bold"), strip.background=element_rect(fill="gray95"),
        panel.border = element_blank(), axis.line=element_line(color="black"),
        legend.position="top", legend.title=element_text(size=20), legend.text=element_text(size=15,color="black"))

BOXPLOT_THEME2 <- BOXPLOT_THEME
BOXPLOT_THEME2$axis.title.x <- element_text(size=20,color="black")

BARPLOT_THEME <- theme_classic() +
  theme(axis.text.x=element_blank(), axis.text.y=element_text(size=15,color="black"),
        axis.title.x=element_text(size=20,color="black"), axis.title.y=element_text(size=20,color="black"),
        panel.border = element_blank(), axis.line=element_line(color="black"),
        legend.position="top", legend.title=element_text(size=20), legend.text=element_text(size=15,color="black"))

