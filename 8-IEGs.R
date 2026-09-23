#Author: Debora Desantis (DD)
#Email: desantis@connect.hku.hk (DD)
#Affilitions: The Hong Kong University (HKU) - The Swire Institute for Marine Sciences (SWIMS)
#Copyright: Copyright 2025, Debora Desantis

#IEG-like Gene Identification & Treatment-Associated IEG Expression Analysis
#Based on: Johnson et al., Nature Communications (2023)

library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
library(lme4)       
library(lmerTest)   
library(broom.mixed)
library(purrr)

#Set seed for reproducibility
set.seed(42)

#Load dataset
Integrated <- readRDS("Integrated_dataset_ctrl-Int_FindAllMarkers_ready_0.8_ClusterNames.rds")

#IEG list: fosb, egr1, npas4
canonical_IEGs <- c("g3812", "g14916", "g21363")

#set params
min_IEG_pos    <- 3
logfc_thresh   <- 0
pval_thresh    <- 0.05
majority_frac  <- 0.5
cluster_sizes <- table(Integrated$CellType)
min_pct_thresh <- 1 / min(cluster_sizes)

#Set default assay
DefaultAssay(Integrated) <- "SCT"

clusters <- levels(factor(Integrated$CellType))
cat("Number of clusters:", length(clusters), "\n")

missing_IEGs <- setdiff(canonical_IEGs, rownames(Integrated))
if (length(missing_IEGs) > 0) {
  warning("The following IEGs were NOT found in the dataset: ",
          paste(missing_IEGs, collapse = ", "))
  canonical_IEGs <- intersect(canonical_IEGs, rownames(Integrated))
}

###---Identify IEG co-expressed genes per cluster per each canonical IEG----
get_expr <- function(obj, gene) {
  GetAssayData(obj, assay = "SCT", layer = "data")[gene, , drop = TRUE]
}

ieg_pval_list <- setNames(vector("list", length(canonical_IEGs)), canonical_IEGs)

for (ieg in canonical_IEGs) {
  
  cat("Processing IEG:", ieg, "\n")

  ieg_expr  <- get_expr(Integrated, ieg)
  all_genes <- rownames(Integrated)
  
  pval_matrix <- matrix(1, nrow = length(all_genes), ncol = length(clusters),
                        dimnames = list(all_genes, clusters))
  lfc_matrix  <- matrix(0, nrow = length(all_genes), ncol = length(clusters),
                        dimnames = list(all_genes, clusters))
  
  clusters_used <- c()
  
  for (cl in clusters) {
    
    cells_in_cluster <- WhichCells(Integrated, expression = CellType == cl)
    ieg_cl           <- ieg_expr[cells_in_cluster]
    
    IEG_pos_cells <- names(ieg_cl[ieg_cl > 0])
    IEG_neg_cells <- names(ieg_cl[ieg_cl == 0])
    
    if (length(IEG_pos_cells) < min_IEG_pos) {
      cat("  Cluster", cl, "skipped: only", length(IEG_pos_cells),
          "IEG+ nuclei\n")
      next
    }
    
    clusters_used <- c(clusters_used, cl)
    cat("  Cluster:", cl,
        "| IEG+:", length(IEG_pos_cells),
        "| IEG-:", length(IEG_neg_cells), "\n")
    
    sub_obj <- subset(Integrated, cells = c(IEG_pos_cells, IEG_neg_cells))
    sub_obj$IEG_status <- ifelse(
      colnames(sub_obj) %in% IEG_pos_cells, "IEG_pos", "IEG_neg")
    Idents(sub_obj) <- "IEG_status"
    
    markers <- tryCatch(
      FindMarkers(
        sub_obj,
        ident.1         = "IEG_pos",
        ident.2         = "IEG_neg",
        logfc.threshold = logfc_thresh,
        min.pct         = min_pct_thresh,
        test.use        = "wilcox",
        assay           = "SCT",
        recorrect_umi   = FALSE),
      error = function(e) {
        warning("FindMarkers failed for IEG: ", ieg,
                " | Cluster: ", cl, " | Error: ", e$message)
        return(NULL)
      }
    )
    
    if (!is.null(markers) && nrow(markers) > 0) {
      tested_genes <- rownames(markers)
      pval_matrix[tested_genes, cl] <- markers$p_val
      lfc_matrix[tested_genes,  cl] <- markers$avg_log2FC
    }
    
  } # end cluster loop

###---Identify genes significantly co-expressed with the IEG----
n_clusters_used <- length(clusters_used)
  
  if (n_clusters_used == 0) {
    warning("No eligible clusters for IEG: ", ieg)
    ieg_pval_list[[ieg]] <- data.frame(gene = character(),
                                       stringsAsFactors = FALSE)
    next
  }
  
  pval_sub <- pval_matrix[, clusters_used, drop = FALSE]
  lfc_sub  <- lfc_matrix[,  clusters_used, drop = FALSE]
  
  n_tested      <- rowSums(pval_sub < 1)
  sig_up_matrix <- (pval_sub < pval_thresh) & (lfc_sub > 0)
  n_sig_up      <- rowSums(sig_up_matrix)
  
  crit_a <- n_tested >= ceiling(majority_frac * n_clusters_used)
  crit_b <- ifelse(n_tested > 0, (n_sig_up / n_tested) > majority_frac, FALSE)
  
  co_expressed_genes <- setdiff(all_genes[crit_a & crit_b], canonical_IEGs)
  
  cat("\n  Genes significantly co-expressed with", ieg, ":",
      length(co_expressed_genes), "\n")
  
  ieg_pval_list[[ieg]] <- data.frame(
    gene       = co_expressed_genes,
    n_clusters = n_tested[co_expressed_genes],
    n_sig_up   = n_sig_up[co_expressed_genes],
    stringsAsFactors = FALSE)
  
} # end IEG loop

###---IEG-like genes when co-expressed with all the three canonical IEGs
co_expr_sets   <- lapply(ieg_pval_list, function(df) df$gene)
IEG_like_genes <- Reduce(intersect, co_expr_sets)

###---IEG Score Analysis with GLMM-Negative Binomial----
IEG_like_genes <- intersect(IEG_like_genes, rownames(Integrated))

#Assign IEG score
counts_matrix <- GetAssayData(Integrated, assay = "RNA", layer = "counts")[IEG_like_genes, ]
ieg_score <- colSums(counts_matrix > 0)

#Add score and log-transformed total features to metadata
Integrated$IEG_score <- ieg_score
Integrated$log_nFeature <- log(Integrated$nFeature_RNA)
Integrated$Treatment <- factor(Integrated$Treatment, levels = c("ctrl", "Int"))
Integrated$sample    <- factor(Integrated$sample)
Integrated$CellType  <- factor(Integrated$CellType)

#Analyze differences per cluster using GLMM (glmer.nb)
results_list <- list()

clusters <- levels(as.factor(Integrated$CellType))

for (cluster in clusters) {
  message("Processing cluster: ", cluster)
  
  #Subset data for the specific cluster
  cluster_data <- subset(Integrated, subset = CellType == cluster)@meta.data
  
  #Fit the negative binomial generalized linear mixed-effects model
  tryCatch({
    model <- glmer.nb(IEG_score ~ Treatment + (1 | sample) + offset(log_nFeature),
                      data = cluster_data,
                      control = glmerControl(optimizer = "bobyqa")) #Optimizer
    
    #Extract coefficients
    summary_mod <- summary(model)$coefficients
    
    if ("TreatmentInt" %in% rownames(summary_mod)) {
      res <- as.data.frame(t(summary_mod["TreatmentInt", ]))
      res$Cluster <- cluster
      results_list[[cluster]] <- res
    }
  }, error = function(e) {
    message("Model failed for cluster ", cluster, ": ", e$message)
  })
}

#Combine results and apply FDR
final_results <- bind_rows(results_list)
final_results$q_val <- p.adjust(final_results$`Pr(>|z|)`, method = "fdr")

#Filter for significance (q < 0.05)
significant_clusters <- final_results %>% filter(q_val < 0.05)