#Author: Debora Desantis (DD)
#Email: desantis@connect.hku.hk (DD)
#Affilitions: The Hong Kong University (HKU) - The Swire Institute for Marine Sciences (SWIMS)
#Copyright: Copyright 2025, Debora Desantis

###hdWGCNA analysis###

library(Seurat)
library(tidyverse)
library(cowplot)
library(patchwork)
library(WGCNA)
library(hdWGCNA)
library(presto)
library(ggplot2)

#Cowplot theme for ggplot
theme_set(theme_cowplot())

#Set seed for reproducibility
set.seed(42)

#Enable multithreading
enableWGCNAThreads(nThreads = 2)

#Load dataset
Integrated <- readRDS("Integrated_dataset_ctrl-Int_FindAllMarkers_ready_0.8_ClusterNames.rds")

#Setup Seurat object for WGCNA
Integrated <- SetupForWGCNA(Integrated,
                            gene_select = "fraction",
                            fraction = 0.05,
                            assay = "SCT",
                            wgcna_name = "Coexpression_Analysis")

#Construct Metacells
Integrated <- MetacellsByGroups(
  seurat_obj = Integrated,
  group.by = c("CellType", "sample"), 
  reduction = "integrated.rpca", 
  k = 25,          
  max_shared = 10, 
  min_cells = 50,
  ident.group = 'CellType')

#Normalize Metacells
Integrated <- NormalizeMetacells(Integrated)

#Set the expression data 
Integrated <- SetDatExpr(Integrated,
                         group_name = "GABA-9", 
                         group.by = "CellType",
                         assay = 'SCT', 
                         layer = 'data')
#Soft powers
Integrated <- TestSoftPowers(Integrated, networkType = 'signed')
plot_list <- PlotSoftPowers(Integrated)
wrap_plots(plot_list, ncol=2)

#Set the soft power
soft_power <- 8

#Construct the network
Integrated <- ConstructNetwork(Integrated,
                               soft_power = soft_power,
                               setDatExpr = FALSE,
                               overwrite_tom = TRUE,
                               tom_name = "GABA-9")
# Visualize the dendrogram
PlotDendrogram(Integrated, main='GABA-9 clusters hdWGCNA Dendrogram')

###---Relate the modules to my Treatment----
#Compute Module Eigengenes (MEs)
Integrated <- ModuleEigengenes(Integrated) #group.by.vars = "sample"

#Extract MEs and merge with metadata for plotting/statistics
#ME with row as cells and columns as modules
MEs <- GetMEs(Integrated, harmonized=FALSE)
metadata <- GetMetacellObject(Integrated)@meta.data
Integrated@meta.data <- cbind(Integrated@meta.data, MEs)

#Compute eigengene-based connectivity (kME)
#Computes pairwise correlations between genes and module eigengenes
Integrated <- ModuleConnectivity(Integrated,
                                 group.by = 'CellType', 
                                 group_name = "GABA-9")

#Rename the modules
Integrated <- ResetModuleNames(Integrated, new_name = "GABA-9-M")

#Plot genes ranked by kME for each module
PlotKMEs(Integrated, ncol=3)

#Get the module assignment table:
modules <- GetModules(Integrated) %>% subset(module != 'grey')

#Get hub genes
hub_df <- GetHubGenes(Integrated, n_hubs = 20)

plot_umap <- ModuleFeaturePlot(Integrated, features = 'hMEs', reduction = "umap.rpca", point_size = 0.1, order = TRUE)
wrap_plots(plot_umap, ncol = 3)

###---DME analysis----
#Identify cell barcodes for each experimental group
group1 <- colnames(Integrated)[Integrated@meta.data$Treatment == "Int"]
group2 <- colnames(Integrated)[Integrated@meta.data$Treatment == "ctrl"]

#Run the Differential Module Eigengene (DME) test
DMEs_df <- FindDMEs(Integrated,
                    barcodes1 = group1,
                    barcodes2 = group2,
                    test.use = "wilcox", 
                    min.pct = 0.25,
                    wgcna_name = "Coexpression_Analysis")

#Generate a volcano plot for the modules
PlotDMEsVolcano(Integrated,
                DMEs_df,
                wgcna_name = "Coexpression_Analysis",
                plot_labels = TRUE)

###Repeat for Glut-12 and GABA-1/GABA-18 combined clusters###

#############################################################
###---scGSEA on GABA-9-M5 module ----
modules <- GetModules(Integrated) %>% subset(module != 'grey')

#rank all genes in this kME
cur_mod <- 'GABA-9-M5'
cur_genes <- modules[,(c('gene_name', 'module', paste0('kME_', cur_mod)))]
ranks <- cur_genes$kME; names(ranks) <- cur_genes$gene_name
ranks <- ranks[order(ranks)]

#rank by GABA-9-M5 genes only by kME
cur_mod <- 'GABA-9-M5'
modules <- GetModules(Integrated) %>% subset(module == cur_mod)
cur_genes <- modules[,(c('gene_name', 'module', paste0('kME_', cur_mod)))]
ranks <- cur_genes$kME; names(ranks) <- cur_genes$gene_name
ranks <- ranks[order(ranks)]

###----scGSEA on hdWGCNA modules----
#Run fgsea
res.gsea <- fgsea(pathways = custom_gene_sets, 
                  stats    = ranks,
                  scoreType = 'std',
                  minSize = 25,
                  maxSize = 500,
                  nproc = 1) 

#res.gsea[is.na(res.gsea)] <- "-"

##Save the results 
data.table::fwrite(res.gsea, file = "gsea_results_module5.tsv", sep = "\t")

sum(res.gsea[, padj < 0.05])
sum(res.gsea[, pval < 0.05])

###Repeat for GABA-9-M1 module###