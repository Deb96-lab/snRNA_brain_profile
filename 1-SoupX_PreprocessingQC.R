#Author: Debora Desantis (DD)
#Email: desantis@connect.hku.hk (DD)
#Affilitions: The Hong Kong University (HKU) - The Swire Institute for Marine Sciences (SWIMS)
#Copyright: Copyright 2025, Debora Desantis

###Data preprocessing, SoupX###

library(Seurat)
library(ggplot2)
library(dplyr)
library(SoupX)
library(sctransform)
library(glmGamPoi)
library(cowplot)

#Set seed for reproducibility
set.seed(42)

###----SoupX to remove ambient RNA from snRNA data. Done per each sample----
###----SoupX Auto
sc = load10X("LD21-ctrl/outs")
sc = autoEstCont(sc)
out = adjustCounts(sc)

#Save results to import with Seurat
#DropletUtils:::write10xCounts("LD22-ctrl/LD22ctrl_SoupX_out", out)

###----Pre-process QC per individual sample----
LD21ctrl_seurat <- CreateSeuratObject(counts = out, project = "LD21ctrl", 
                                      min.cells = 3)
LD21ctrl_seurat

##Add the %MT to the metadata
LD21ctrl_seurat[["percent.mt"]] <- PercentageFeatureSet(LD21ctrl_seurat, 
                                                 features = c("COXI", "COXII", "COXIII", "ND1", "ND2", 
                                                              "ND3", "ND4", "ND5", "ND6", "ND4L",
                                                              "ATPase8", "ATPase6", "Cytb"))

###----Generate QC plots to inspect data---- 
plot1 <- VlnPlot(LD21ctrl_seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
                 group.by = "orig.ident", ncol = 3)
plot1 

#plot the number of counts against the number of features (genes)
plot2 <- FeatureScatter(LD21ctrl_seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot2

#Plot the number of counts against the number of MT genes
plotMT <- FeatureScatter(LD21ctrl_seurat, feature1 = "nCount_RNA", feature2 = "percent.mt")
plotMT

###----Filter dataset----
LD21ctrl_seurat <- subset(LD21ctrl_seurat, subset = nFeature_RNA > 250 & nFeature_RNA < 5000 &
                           nCount_RNA < 18000 & percent.mt < 3)
LD21ctrl_seurat

#Plot again to check the data distribution
plot3 <- VlnPlot(LD21ctrl_seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
                 ncol = 3, pt.size = 0.01, group.by = "orig.ident")
plot3

###Repeat for each individual sample and move to scDblFinder###