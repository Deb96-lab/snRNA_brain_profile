#Author: Debora Desantis (DD)
#Email: desantis@connect.hku.hk (DD)
#Affilitions: The Hong Kong University (HKU) - The Swire Institute for Marine Sciences (SWIMS)
#Copyright: Copyright 2025, Debora Desantis

###Re-cluster RGc cell clusters###

#Load required libraries
library(readxl)
library(stringr)
library(tibble)
library(tidyr)
library(dplyr)
library(Seurat)
library(ggplot2)
library(patchwork)
library(lme4) 
library(tidyverse)

#Set seed for reproducibility
set.seed(42)

#Load dataset
Integrated <- readRDS("Integrated_dataset_ctrl-Int_FindAllMarkers_ready_0.8_ClusterNames.rds")

#Add metadata column for CellType
Integrated$CellType <- Idents(Integrated)

Integrated$CellType <- gsub("^RGc-.*", "RGc", Integrated$CellType)

#Subset 
Idents(Integrated) <- "CellType"
RGc_sub <- subset(Integrated, idents = "RGc")

#Split layers
RGc_sub[["RNA"]] <- split(RGc_sub[["RNA"]], f = RGc_sub$orig.ident)
Layers(RGc_sub[["RNA"]])
rm(Integrated)

###---Normalize with SCTransform and identify variable features per each dataset independently----
RGc_sub <- SCTransform(RGc_sub, method = "glmGamPoi", 
                          vars.to.regress = c("percent.mt"), verbose = TRUE) 

#Run PCA
RGc_sub <- RunPCA(RGc_sub)

#Join layers 
RGc_sub[["RNA"]] <- JoinLayers(RGc_sub[["RNA"]])
Layers(RGc_sub[["RNA"]])

#Cluster Cells to find the one with highest ribo %
RGc_sub <- FindNeighbors(RGc_sub, reduction = "integrated.rpca", assay = "SCT",
                            dims = 1:30, verbose = TRUE)
RGc_sub <- FindClusters(RGc_sub, resolution = 0.5, cluster.name = "RGc_umap", 
                           RGc_subverbose = TRUE) 
RGc_sub <- RunUMAP(RGc_sub, reduction = "integrated.rpca", reduction.name = "UMAP",
                      dims = 1:30, verbose = TRUE)

#Rename Clusters
RGc_sub <- RenameIdents(object = RGc_sub, `0` = "RGc-0")
RGc_sub <- RenameIdents(object = RGc_sub, `1` = "RGc-1")
RGc_sub <- RenameIdents(object = RGc_sub, `2` = "RGc-2")
RGc_sub <- RenameIdents(object = RGc_sub, `3` = "RGc-3")
RGc_sub <- RenameIdents(object = RGc_sub, `4` = "RGc-4")
RGc_sub <- RenameIdents(object = RGc_sub, `5` = "RGc-5")
RGc_sub <- RenameIdents(object = RGc_sub, `6` = "RGc-6")

metadata_df <- RGc_sub@meta.data
RGc_sub$CellType <- Idents(RGc_sub)

################################################################################
###---Monocle3 Trajectory analysis all clusters----
library(monocle3)
library(SeuratWrappers)

cds <- as.cell_data_set(RGc_sub)
cds <- cluster_cells(cds, reduction_method = "UMAP")

p1 <- plot_cells(cds, color_cells_by = "partition", show_trajectory_graph = FALSE)

###---Trajectory analysis on partition 1----
#Identify which cells belong to partition 1
cds <- learn_graph(cds, use_partition = TRUE)

p2 <- plot_cells(cds, color_cells_by = "CellType",
                 label_groups_by_cluster = FALSE, 
                 label_leaves = FALSE, 
                 label_branch_points = FALSE,
                 show_trajectory_graph = FALSE)

#cds <- order_cells(cds)

#Find the root (biological a priopri knowledge)
get_earliest_principal_node <- function(cds,
                                        root_group = "RGc-1",
                                        group_col = "CellType") {
  cell_ids <- which(colData(cds)[[group_col]] == root_group)
  
  #Map each cell to its closest vertex on the principal graph
  closest_vertex <- cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  
  #Principal graph nodes
  root_pr_nodes <- igraph::V(principal_graph(cds)[["UMAP"]])$name
  root_pr_nodes[
    as.numeric(names(which.max(table(closest_vertex[cell_ids, ]))))
  ]
}

#Order cells 
cds <- order_cells(cds,
                   root_pr_nodes = get_earliest_principal_node(cds, "RGc-1", "CellType"))