#Author: Debora Desantis (DD)
#Email: desantis@connect.hku.hk (DD)
#Affilitions: The Hong Kong University (HKU) - The Swire Institute for Marine Sciences (SWIMS)
#Copyright: Copyright 2025, Debora Desantis

#Load required libraries
library(Seurat)
library(tidyverse)
library(CellChat)
library(patchwork)
options(stringsAsFactors = FALSE)

#Set seed for reproducibility
set.seed(42)

#Load dataset
Integrated <- readRDS("Integrated_dataset_ctrl-Int_FindAllMarkers_ready_0.8_ClusterNames.rds")

###Merge RGc sub info with Integrated obj----
metadata_df <- RGc_sub@meta.data
RGc_sub$CellType <- Idents(RGc_sub)
subcluster_labels <- RGc_sub$CellType

#Add these labels as a new metadata column in the main object
Integrated <- AddMetaData(
  object = Integrated,
  metadata = subcluster_labels,
  col.name = "subclusters")

#Convert the new metadata column to a character vector
current_labels <- as.character(Integrated$subclusters)

#Fix NA values
na_indices <- which(is.na(current_labels))
current_labels[na_indices] <- as.character(Integrated$CellType[na_indices])
Integrated$combined_clusters <- factor(current_labels)

Idents(Integrated) <- "combined_clusters"

###Add mouse annotation----
ortho_map <- read.csv("orthologs_with_genes_mouse.csv")

#Extract raw counts from your Seurat object
exp_matrix <- GetAssayData(Integrated, assay = "RNA", layer = "counts")

#Filter and Rename Genes based on Orthologs  to keep only genes present in mapping file
genes_to_keep <- intersect(rownames(exp_matrix), ortho_map$Ldim_ID)
exp_matrix <- exp_matrix[genes_to_keep, ]

#Create a named vector for remapping
mapping_vector <- ortho_map$Gene_Name
names(mapping_vector) <- ortho_map$Ldim_ID

#Update row names to Mouse Symbols
new_rownames <- mapping_vector[rownames(exp_matrix)]
exp_matrix <- exp_matrix[!is.na(new_rownames), ]  
rownames(exp_matrix) <- new_rownames[!is.na(new_rownames)]

#Handle duplicates before creating Seurat obj
if (any(duplicated(rownames(exp_matrix)))) {
  message("Duplicate gene symbols detected. Aggregating counts...")
  exp_matrix <- rowsum(
    as.matrix(exp_matrix),
    group = rownames(exp_matrix))
}

###---Create a clean Seurat object for CellChat----
seurat_object <- CreateSeuratObject(counts = exp_matrix)
seurat_object@meta.data <- Integrated@meta.data[colnames(seurat_object), ]

#Split layers
seurat_object[["RNA"]] <- split(seurat_object[["RNA"]], f = seurat_object$orig.ident)
Layers(seurat_object[["RNA"]])
rm(Integrated, exp_matrix, RGc_sub)
gc()

#Run SCTransform pipeline
###---Normalize with SCTransform and identify variable features per each dataset independently
seurat_object <- SCTransform(seurat_object, method = "glmGamPoi", 
                             vars.to.regress = c("percent.mt"), verbose = TRUE) 

#Run PCA
seurat_object <- RunPCA(seurat_object, verbose = TRUE)
gc()

#Join layers 
seurat_object[["RNA"]] <- JoinLayers(seurat_object[["RNA"]])
Layers(seurat_object[["RNA"]])

###----Start CellChat Workflow----
#Subset 
Idents(seurat_object) <- "combined_clusters"

data.input <- seurat_object[["SCT"]]@data
labels <- Idents(seurat_object)
meta <- data.frame(labels = labels, row.names = names(labels))

seurat_object$samples <- seurat_object$orig.ident
cellChat <- createCellChat(object = seurat_object, group.by = "ident", assay = "SCT")

CellChatDB <- CellChatDB.mouse 
#showDatabaseCategory(CellChatDB)
CellChatDB.use <- CellChatDB 

#Assign the database to the cellChat object
cellChat@DB <- CellChatDB.use
cellChat <- subsetData(cellChat)
cellChat <- identifyOverExpressedGenes(cellChat)
cellChat <- identifyOverExpressedInteractions(cellChat)

ptm = Sys.time()
cellChat <- computeCommunProb(cellChat, type = "triMean", population.size = FALSE)

cellChat <- filterCommunication(cellChat, min.cells = 10)

cellChat <- computeCommunProbPathway(cellChat)

cellChat <- aggregateNet(cellChat)

#Extracting the inferred communications at the level of ligands/receptors
df.net <- subsetCommunication(cellChat)

#Extracting the inferred communications at the level of signaling pathways
df.netP <- subsetCommunication(cellChat, slot.name = "netP")