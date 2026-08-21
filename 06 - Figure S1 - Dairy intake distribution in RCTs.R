# Date : 16 septembre 2025
# Auteur : Jacob Lessard-Lord

# Projet PLC

# Figure S1 - Milk and cheese intake distribution in fully controlled trials

# Set working directory ----
setwd("C:/Users/jackl/OneDrive/Bureau/École/Postdoctorat/Chercheur/Jean-Philippe/PLC/Processing_PLC/Untargeted/Final_Processing/Data")

# Chargement des packages requis ----
library(tidyverse)
library(ggpubr)

# Importation des données métabolomiques ----
load("Metabo/Untargeted - ALL.rda")

# Importation des données cliniques ----
load("FANI/Data clinique - RCT.rda")

# Mise en forme des données ----

# Selection des données post-dietes RCT
data_ech <- subset(data_ech, subset = Traitement %in% c("Lait", "Fromage", "Beurre", "Controle"),
                   select = c("Projet", "Sujet", "Projet_Sujet", "Traitement"))

# Merge avec les données cliniques

## Lait 
data_rct_lait <- merge(x = data_rct_lait %>% mutate(Traitement = factor(str_replace(as.character(Traitement), "Contrôle", "Controle"))),
                       data_ech, by = c("Projet", "Sujet", "Traitement")) %>%
  filter(Projet_Sujet != "LAIT_16")

## Fromage
data_rct_fromage <- merge(x = data_rct_fromage %>% mutate(Traitement = factor(str_replace(as.character(Traitement), "Contrôle", "Controle"))),
                          data_ech, by = c("Projet", "Sujet", "Traitement"))

## Beurre
data_rct_beurre <- merge(x = data_rct_beurre %>% mutate(Traitement = factor(str_replace(as.character(Traitement), "Contrôle", "Controle"))),
                          data_ech, by = c("Projet", "Sujet", "Traitement"))

# Mise en forme pour graphique
data_graph <- rbind(data_rct_lait %>% filter(Traitement == "Lait") %>%
                      select("Milk_intake") %>%
                      rename("Dairy_intake" = "Milk_intake") %>%
                      mutate(Food = "Milk"),
                    data_rct_fromage %>%
                      filter(Traitement == "Fromage") %>%
                      select("Cheese_intake") %>%
                      rename("Dairy_intake" = "Cheese_intake") %>%
                      mutate(Food = "Cheese"),
                    data_rct_beurre %>%
                      filter(Traitement == "Beurre") %>%
                      select("Butter_intake") %>%
                      rename("Dairy_intake" = "Butter_intake") %>%
                      mutate(Food = "Butter")) %>%
  mutate(Food = factor(Food, levels = c("Milk", "Cheese", "Butter")))

# Obtenir statistiques descriptives de la distribution
data_graph %>% group_by(Food) %>% get_summary_stats(Dairy_intake, type = "full")

# Graphique ----
setwd("C:/Users/jackl/OneDrive/Bureau/École/Postdoctorat/Chercheur/Jean-Philippe/PLC/Article/Biomarqueurs/Figures")

tiff(filename = "Figure S1 - Distribution Dairy intake - RCT.tiff", width = 15, height = 22.5, units = "cm", res = 600, compression = "lzw")

gghistogram(data = data_graph, x = "Dairy_intake",
            xlab = "Daily intake (g)", ylab = "Frequency",
            facet.by = "Food", nrow = 3, scales = "free_x", fill = "Food", palette = c("#425ae1", "#fd533d", "#FFD166"),
            panel.labs.font = list(size = 11, face = "bold", color = "white"), panel.labs.background = list(fill = "black")) +
  scale_y_continuous(breaks = c(0, 2, 4, 6, 8, 10)) +
  theme(legend.position = "none")

dev.off()
