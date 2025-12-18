# Date : 21 novembre 2025
# Auteur : Jacob Lessard-Lord

# Projet PLC

# Figure 2 - Sélection des features par XGBoost 

# Set working directory ----
setwd("C:/Users/jackl/OneDrive/Bureau/École/Postdoctorat/Chercheur/Jean-Philippe/PLC/Processing_PLC/Untargeted/Final_Processing/Data")

# Chargement des packages requis ----
library(tidyverse)
library(ggpubr)

# Importation des données ----
load("Metabo/XGB - Selection variables.rda")

# Mise en forme des données ----
results_xgb_select <- rbind(results_xgb_select_lait %>% mutate(Classification = "Milk vs No dairy"),
                            results_xgb_select_fromage %>% mutate(Classification = "Cheese vs No dairy"),
                            results_xgb_select_beurre %>% mutate(Classification = "Butter vs No dairy"),
                            results_xgb_select_lait_fromage %>% mutate(Classification = "Milk vs Cheese"),
                            results_xgb_select_fromage_beurre %>% mutate(Classification = "Cheese vs Butter")) %>%
  pivot_longer(cols = c("AUC", "Sensibilite", "Specificite"), names_to = "Metrics", values_to = "Valeur") %>%
  mutate(Metrics = case_match(Metrics, "AUC" ~ "AUROC", "Sensibilite" ~ "Sensibility", "Specificite" ~ "Specificity"),
         Classification = factor(Classification, levels = c("Milk vs No dairy", "Cheese vs No dairy", "Butter vs No dairy", "Milk vs Cheese", "Cheese vs Butter")))

# Graphique ----

setwd("C:/Users/jackl/OneDrive/Bureau/École/Postdoctorat/Chercheur/Jean-Philippe/PLC/Article/Biomarqueurs/Figures")

tiff(filename = "Figure 2 - Perfomance XGBoost.tiff", width = 30, height = 20, units = "cm", res = 600, compression = "lzw")

ggline(data = results_xgb_select, x = "Split", y = "Valeur", add = "mean_ci",
       facet.by = c("Metrics", "Classification"), scales = "free_x",
       panel.labs.font = list(size = 11, face = "bold", color = "white"), panel.labs.background = list(fill = "black"),
       xlab = "Iteration") +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "red") +
  theme(axis.title.y = element_blank())

dev.off()