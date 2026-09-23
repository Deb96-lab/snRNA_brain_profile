#Author: Debora Desantis (DD)
#Email: desantis@connect.hku.hk (DD)
#Affilitions: The Hong Kong University (HKU) - The Swire Institute for Marine Sciences (SWIMS)
#Copyright: Copyright 2025, Debora Desantis

###FindAllMarkers###

library(Seurat)
library(ggplot2)
library(dplyr)
library(sctransform)
library(glmGamPoi)
library(presto)

parallel::mclapply

#Set seed for reproducibility
set.seed(42)

###----Cluster annotation----
#Load data
Integrated <- readRDS("Integrated_dataset_ctrl-Int_clusterInfo_rpca_0.8.rds")

Layers(Integrated[["RNA"]])

##Prep dataset for FindAllMarkers function
Integrated <- PrepSCTFindMarkers(Integrated, assay = "SCT")

##Find markers
Int_markers <- FindAllMarkers(Integrated, only.pos = TRUE, 
                              min.pct = 0.25, logfc.threshold = 0.25, test.use = "wilcox")
#Save dataset
write.csv(Int_markers, file = "outs/Integrated_ctrl-Int_markers_rpca-0.8.csv")

#Select the top 5 markers and save the output
Top5 <- Int_markers %>%
  group_by(cluster) %>%
  slice_max(n = 5, order_by = avg_log2FC)

write.csv(Top5, file = "outs/Integrated_ctrl-Int_markers_rpca_top5.csv")

#Plot Top5 markers from FindAllMarkers
d1 <-DotPlot(Integrated, features = Top5, cols = c("lightgrey", "blue"),
             col.min = 0, dot.scale = 6) +  
  theme_bw(base_size = 8, base_family = plot_font) +
  scale_x_discrete() +  
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 6, color = "black"),
    axis.text.y = element_text(angle = 0, hjust = 1, size = 6, color = "black"),
    axis.title.x = element_text(size = 8, face = "bold", color = "black", margin = margin(t = 8)),
    axis.title.y = element_text(size = 8, face = "bold", color = "black",  margin = margin(r = 8)),
    legend.title = element_text(size = 8, angle = 0),
    legend.text  = element_text(size = 6, color = "black"),
    panel.grid   = element_blank(),
    axis.ticks   = element_blank()) +
  guides(fill = guide_colourbar(
    title.position = "top",
    barheight = unit(1.5, "cm"),
    barwidth  = unit(0.3, "cm"))) + 
  coord_flip()

###---Rename clusters----
Integrated <- RenameIdents(object = Integrated, `0` = "Glut-0")
Integrated <- RenameIdents(object = Integrated, `1` = "GABA-1")
Integrated <- RenameIdents(object = Integrated, `2` = "Glut-2")
Integrated <- RenameIdents(object = Integrated, `3` = "Glut-3")
Integrated <- RenameIdents(object = Integrated, `4` = "Pallium-4")
Integrated <- RenameIdents(object = Integrated, `5` = "Glut-5")
Integrated <- RenameIdents(object = Integrated, `6` = "Glut-6")
Integrated <- RenameIdents(object = Integrated, `7` = "Glut-7")
Integrated <- RenameIdents(object = Integrated, `8` = "Glut-8")
Integrated <- RenameIdents(object = Integrated, `9` = "GABA-9")
Integrated <- RenameIdents(object = Integrated, `10` = "RGc-10")
Integrated <- RenameIdents(object = Integrated, `11` = "Glut-11")
Integrated <- RenameIdents(object = Integrated, `12` = "Glut-12")
Integrated <- RenameIdents(object = Integrated, `13` = "GABA-13")
Integrated <- RenameIdents(object = Integrated, `14` = "RGc-14")
Integrated <- RenameIdents(object = Integrated, `15` = "Glut-15")
Integrated <- RenameIdents(object = Integrated, `16` = "Glut-16")
Integrated <- RenameIdents(object = Integrated, `17` = "Oligo")
Integrated <- RenameIdents(object = Integrated, `18` = "GABA-18")
Integrated <- RenameIdents(object = Integrated, `19` = "Glut-19")
Integrated <- RenameIdents(object = Integrated, `20` = "Glut-20")
Integrated <- RenameIdents(object = Integrated, `21` = "OPCs")
Integrated <- RenameIdents(object = Integrated, `22` = "Perithelial")

#Add metadata column for CellType
Integrated$CellType <- Idents(Integrated)

saveRDS(Integrated, "Integrated_dataset_ctrl-Int_FindAllMarkers_ready_0.8_ClusterNames.rds")