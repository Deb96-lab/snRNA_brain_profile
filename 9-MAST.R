#Author: Debora Desantis (DD)
#Email: desantis@connect.hku.hk (DD)
#Affilitions: The Hong Kong University (HKU) - The Swire Institute for Marine Sciences (SWIMS)
#Copyright: Copyright 2025, Debora Desantis

###MAST Differential Expression Analysis###

library(Seurat)
library(SingleCellExperiment)
library(MAST)
library(dplyr)
library(data.table)
library(future)
library(future.apply)
library(dplyr)
library(purrr)
library(tidyr)
library(ggplot2)

#Set seed for reproducibility
set.seed(42)

#Load dataset
Integrated <- readRDS("Integrated_dataset_ctrl-Int_FindAllMarkers_ready_0.8_ClusterNames.rds")

#Create output directories
dir.create("output/raw_data_for_MAST", recursive = TRUE, showWarnings = FALSE)
dir.create("output/DE_MAST",           recursive = TRUE, showWarnings = FALSE)

#Set factor levels
Integrated$Treatment <- factor(Integrated$Treatment, levels = c("ctrl", "Int"))
Integrated$sample    <- factor(Integrated$sample,
                               levels = c("LDctrl21", "LD22ctrl", "LD31ctrl",
                                          "LD24Int",  "LD26Int",  "LD28Int"))

#Set default assay
DefaultAssay(Integrated) <- "RNA"

###---Split object by CellType and save SCE files----
Idents(Integrated) <- "CellType"
cell_types <- as.character(unique(Idents(Integrated)))

for (ct in cell_types) {
  sub     <- subset(Integrated, idents = ct)
  sce     <- Seurat::as.SingleCellExperiment(sub, assay = "RNA")
  ct_safe <- gsub("[^A-Za-z0-9_.-]", "_", ct)   # sanitise for filename
  saveRDS(sce, paste0("output/raw_data_for_MAST/", ct_safe, "_data.rds"))
  message("  Saved SCE for: ", ct)
}

rm(Integrated, sub, sce)
gc()

###---MAST wrapper function----
calc_mast_DE <- function(ct_safe) {
  sce <- readRDS(paste0("output/raw_data_for_MAST/", ct_safe, "_data.rds"))
  
  assayNames(sce) 
  logcounts(sce) <- log2(counts(sce) + 1) #log2 transform the raw data
  
  #remove undetected genes
  sce <- sce[rowSums(counts(sce) > 0) > 0, ]
  dim(sce)
  
  #remove lowly expressed genes expressed < 10% of cells
  sce <- sce[rowSums(counts(sce) > 1) >= 0.1*ncol(sce), ]
  dim(sce)
  
  #filter out non-protein coding and artefact-prone genes (MT genes)
  nonsense <- c(rownames(sce[grep("^COX", rownames(sce)),]),
              rownames(sce[grep("^ND", rownames(sce)),]),
              rownames(sce[grep("16S-rRNA", rownames(sce)),]),
              rownames(sce[grep("12S-rRNA", rownames(sce)),]),
              rownames(sce[grep("Cytb", rownames(sce)),]))
  sce <- sce[!rownames(sce) %in% nonsense,]
  dim(sce)
  
  sca <- SceToSingleCellAssay(sce, class = "SingleCellAssay")
  
  #add a column to the data which contains scaled number of genes that are expressed in each cell
  cdr2 <- colSums(assay(sca)>0)
  colData(sca)$ngeneson <- scale(cdr2)
  colData(sca)$mt_reads <- scale(colData(sca)$percent.mt)
  colData(sca)$orig.ident <- factor(colData(sca)$orig.ident)
  colData(sca)$Treatment <- factor(colData(sca)$Treatment)
  
  zlmCond <- zlm(formula = ~ngeneson + mt_reads + Treatment + (1 | orig.ident), 
                 sca=sca, 
                 method='glmer', 
                 ebayes=F,
                 parallel = FALSE,
                 silent = TRUE,
                 strictConvergence=F,
                 fitArgsD=list(nAGQ = 0))
  
  #perform LRT 
  summaryCond <- summary(zlmCond, doLRT="TreatmentInt")
  
  #get the table with log-fold changes and p-values
  summaryDt <- summaryCond$datatable
  result <- merge(summaryDt[contrast=='TreatmentInt' & component=='H',.(primerid, `Pr(>Chisq)`)], # p-values
                  summaryDt[contrast=='TreatmentInt' & component=='logFC', .(primerid, coef, ci.hi, ci.lo)],
                  by='primerid') # logFC coefficients
  
  #MAST uses natural logarithm so we convert the coefficients to log2 base
  result[,coef:=result[,coef]/log(2)]
  
  #Multiple testing correction
  result[,FDR:=p.adjust(`Pr(>Chisq)`, 'fdr')]
  result <- stats::na.omit(as.data.frame(result))
  write.csv(result,
            paste0("output/DE_MAST/", ct_safe, ".csv"),
            row.names = FALSE)
  
  return(result)
}
 
###---Sequential loop with per-cell-type----
cell_types_safe <- gsub("[^A-Za-z0-9_.-]", "_", cell_types)
 
de_results <- setNames(vector("list", length(cell_types_safe)), cell_types_safe)
 
for (ct in cell_types_safe) {
  de_results[[ct]] <- tryCatch(
    calc_mast_DE(ct),
    error = function(e) {
      message("  ERROR in ", ct, ": ", conditionMessage(e))
      NULL
    }
  )
}
 
##Combine all results
clean_results <- Filter(Negate(is.null), de_results)
de_combined <- dplyr::bind_rows(clean_results, .id = "CellType")

#Significant results only
sig_de_combined <- de_results %>%
  keep(~ is.data.frame(.x)) %>%      
  bind_rows(.id = "cell_type") %>%    
  filter(FDR < 0.05)                  