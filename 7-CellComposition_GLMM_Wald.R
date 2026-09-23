#Author: Debora Desantis (DD)
#Email: desantis@connect.hku.hk (DD)
#Affilitions: The Hong Kong University (HKU) - The Swire Institute for Marine Sciences (SWIMS)
#Copyright: Copyright 2025, Debora Desantis

###GLMM with Wald Test for cell composition analysis###

library(lme4)
library(dplyr)
library(purrr)
library(broom.mixed) 
library(Matrix)
library(tibble)
library(ggplot2)
library(patchwork)

#Set seed for reproducibility
set.seed(42)

#Load data
Integrated <- readRDS("Integrated_dataset_ctrl-Int_FindAllMarkers_ready_0.8_ClusterNames.rds")

#Extract metadata and set ctrl the reference 
meta_all <- Integrated@meta.data %>%
  mutate(Treatment = factor(Treatment, levels = c("ctrl", "Int")),
         sample = as.factor(sample))

clusters <- unique(meta_all$CellType)

#Iterate the analysis through each cluster
stat_results <- map_dfr(clusters, function(cl) {
  
  df_cl <- meta_all %>%
    mutate(is_target = ifelse(CellType == cl, 1, 0))
  result <- tryCatch({
    
    #Fit Model
    model_full <- glmer(is_target ~ Treatment + (1|sample), 
                        data = df_cl, 
                        family = binomial,
                        control = glmerControl(optimizer = "bobyqa"))
    
    #Extract the coefficients table (Wald Test results)
    coef_summary <- summary(model_full)$coefficients
    log_OR <- coef_summary["TreatmentInt", "Estimate"]
    p_val_wald <- coef_summary["TreatmentInt", "Pr(>|z|)"]
    
    tibble(
      CellType = cl,
      log_odds_ratio = log_OR,
      p_value = p_val_wald)
    
  }, error = function(e) {
    return(tibble(CellType = cl, p_value = NA, log_odds_ratio = NA))
  })
  
  return(result)
})
#Loop ended

#Multiple Testing Correction
stat_results <- stat_results %>%
  filter(!is.na(p_value)) %>%
  mutate(
    q_value = p.adjust(p_value, method = "BH"),
    significant = q_value < 0.05,
    direction = case_when(
      log_odds_ratio > 0 & significant ~ "Increased in Int",
      log_odds_ratio < 0 & significant ~ "Decreased in Int",
      TRUE ~ "No significant change")) %>%
  arrange(q_value)