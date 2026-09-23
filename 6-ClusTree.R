#Author: Debora Desantis (DD)
#Email: desantis@connect.hku.hk (DD)
#Affilitions: The Hong Kong University (HKU) - The Swire Institute for Marine Sciences (SWIMS)
#Copyright: Copyright 2025, Debora Desantis

###cluster tree reconstruction analysis###
#Unbiased reconstruction of the transcriptional relationships of the different clusters

#Load required libraries
library(tidyr)
library(dplyr)
library(patchwork)
library(Seurat)
library(ape)
library(tidyverse)
library(ggtree)
library(MatrixGenerics)
library(RColorBrewer)

#Set seed for reproducibility
set.seed(42)

#Load dataset
Integrated <- readRDS("Integrated_dataset_ctrl-Int_FindAllMarkers_ready_0.8_ClusterNames.rds")

#Set identity to CellType
Idents(Integrated) <- "CellType"

#Get Average Expression 
avg_matrix <- AverageExpression(Integrated, assays = "SCT", slot = "data", return.seurat = FALSE)

#Identify expressed genes 
expressed_genes <- rownames(avg_matrix$SCT)[which(rowMaxs(as.matrix(avg_matrix$SCT)) > 0.03)]

#Create the final matrix for the tree and subset matrix to expressed gene
clusters <- levels(Integrated)
tree_matrix <- as.matrix(avg_matrix$SCT[expressed_genes, ])

#Helper function for Spearman Distance
calc_spearman_dist <- function(x) {
  as.dist(1 - cor(t(x), method = "spearman"))
}

#Calculate the NJ tree
nj_tree <- nj(calc_spearman_dist(t(tree_matrix)))

#Bootstrapping 
nj_boot <- boot.phylo(phy = nj_tree, 
                      x = t(tree_matrix), 
                      FUN = function(xx) nj(calc_spearman_dist(xx)), 
                      B = 10000) 

#Assign bootstrap values to nodes
nj_tree$node.label <- round(nj_boot / 100)

#Create a summary table for the clusters 
node_metadata <- Integrated@meta.data %>%
  group_by(CellType) %>%
  summarize(
    cell_n = n(), 
    nCount_avg = mean(nCount_RNA), 
    nFeature_avg = mean(nFeature_RNA)) 

###---ggtree plots----
#Initialize ggtree object
gt <- ggtree(nj_tree, layout = "rectangular")

#tree with Bootstrap values
plot_bootstrap <- gt %<+% node_metadata + 
  geom_text2(aes(subset = !isTip, label = label), nudge_x = -0.01, nudge_y = 0.2, size = 3) +
  geom_tiplab(size = 4, colour = "black") +
  ggplot2::xlim(0, max(gt$data$x) + 0.1) + 
  ggtitle("NJ Tree with Bootstrap Support")