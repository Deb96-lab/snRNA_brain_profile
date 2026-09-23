#Author: Debora Desantis (DD)
#Email: desantis@connect.hku.hk (DD)
#Affilitions: The Hong Kong University (HKU) - The Swire Institute for Marine Sciences (SWIMS)
#Copyright: Copyright 2025, Debora Desantis

###Cluster analysis###

library(Seurat)
library(ggplot2)
library(dplyr)
library(sctransform)
library(glmGamPoi)
library(presto)
library(clustree)

parallel::mclapply

#Set seed for reproducibility
set.seed(42)

###----Cluster analysis and visualization of integrated data----
##---Load RPCA Integrated data----
Integrated <- readRDS("Integrated_dataset_ctrl-Int_rpca.rds")

#Apply a standard Seurat workflow
#Join layers 
Integrated[["RNA"]] <- JoinLayers(Integrated[["RNA"]])
Layers(Integrated[["RNA"]])

#Cluster cells RPCA Integration
Integrated <- FindNeighbors(Integrated, reduction = "integrated.rpca", assay = "SCT",
                            dims = 1:30, verbose = TRUE)

###---Run FindCluster at different resolutions
Integrated <- FindClusters(Integrated, resolution = 0.1, cluster.name = "rpca_clusters_0.1", 
                           verbose = TRUE)
Integrated <- FindClusters(Integrated, resolution = 0.2, cluster.name = "rpca_clusters_0.2", 
                           verbose = TRUE)
Integrated <- FindClusters(Integrated, resolution = 0.3, cluster.name = "rpca_clusters_0.3", 
                           verbose = TRUE)
Integrated <- FindClusters(Integrated, resolution = 0.4, cluster.name = "rpca_clusters_0.4", 
                           verbose = TRUE)
Integrated <- FindClusters(Integrated, resolution = 0.5, cluster.name = "rpca_clusters_0.5", 
                           verbose = TRUE)
Integrated <- FindClusters(Integrated, resolution = 0.6, cluster.name = "rpca_clusters_0.6", 
                           verbose = TRUE)
Integrated <- FindClusters(Integrated, resolution = 0.8, cluster.name = "rpca_clusters_0.8", 
                           verbose = TRUE)
Integrated <- FindClusters(Integrated, resolution = 1, cluster.name = "rpca_clusters_1", 
                           verbose = TRUE)
Integrated <- FindClusters(Integrated, resolution = 1.2, cluster.name = "rpca_clusters_1.2", 
                           verbose = TRUE)
Integrated <- FindClusters(Integrated, resolution = 1.3, cluster.name = "rpca_clusters_1.3", 
                           verbose = TRUE)
Integrated <- FindClusters(Integrated, resolution = 1.4, cluster.name = "rpca_clusters_1.4", 
                           verbose = TRUE)
Integrated <- FindClusters(Integrated, resolution = 1.5, cluster.name = "rpca_clusters_1.5", 
                           verbose = TRUE)
Integrated <- FindClusters(Integrated, resolution = 1.6, cluster.name = "rpca_clusters_1.6", 
                           verbose = TRUE)

#Plot the relationship between clustering at different resolution 
c <- clustree(Integrated, prefix = "rpca_clusters_", dims = 1:30) ##res 0.8 is the most appropriate

Integrated <- RunUMAP(Integrated, reduction = "integrated.rpca", reduction.name = "umap.rpca",
                      dims = 1:30, verbose = TRUE)

#Plot UMAP 
p <- DimPlot(Integrated, reduction = "umap.rpca", split.by = "Treatment",
             label = TRUE, repel =TRUE)

#Save integrated dataset with clusters info
saveRDS(Integrated, "Integrated_dataset_ctrl-Int_clusterInfo_rpca_0.8.rds")