# Date : 31 juillet 2025
# Auteur : Jacob Lessard-Lord

# Projet PLC

# ÉValuation de la performance des biomarqueurs en contexte observationnel

# Régression linéaire multiple

# Untargeted 

# Set working directory ----
setwd("C:/Users/jackl/OneDrive/Bureau/École/Postdoctorat/Chercheur/Jean-Philippe/PLC/Processing_PLC/Untargeted/Final_Processing/Data")

# Chargement des packages requis ----
library(tidyverse)
library(ggpubr)
library(readxl)
library(car)
library(olsrr)
library(ComplexHeatmap)
library(circlize)

# Importation des données nutritionnelles ----
load("FFQ/Data FFQ.rda")
load("R24W/Data R24W.rda")

# Importation des données métabolomiques ----
load("Metabo/Untargeted - ALL.rda")

# Importation  des biomarqueurs ----
load("Metabo/Biomarqueurs finaux.rda")

# Importation des scores multi-métabolites ----
load("Metabo/Score multi-metabolite.rda")

# Importation des identifications
identification <- rbind(
  # Lait #
  read_excel("Metabo/ID Biomarqueurs.xlsx", sheet = "Lait") %>% select(c(Chromatography, Polarity, RT, mz, Identification)) %>%
    mutate(ID = paste(Chromatography, Polarity, RT, mz, sep = "_"),
           label = if_else(is.na(Identification) | Identification %in% c("NA",""),
                           ID, Identification)) %>%
    transmute(label, ID) %>%
    filter(ID %in% biomarqueur_lait),
  
  # Beurre #
  read_excel("Metabo/ID Biomarqueurs.xlsx", sheet = "Beurre") %>% select(c(Chromatography, Polarity, RT, mz, Identification)) %>%
    mutate(ID = paste(Chromatography, Polarity, RT, mz, sep = "_"),
           label = if_else(is.na(Identification) | Identification %in% c("NA",""),
                           ID, Identification)) %>%
    transmute(label, ID) %>%
    filter(ID %in% biomarqueur_beurre),
  
  # Dairy fat #
  read_excel("Metabo/ID Biomarqueurs.xlsx", sheet = "Dairy_fat") %>% select(c(Chromatography, Polarity, RT, mz, Identification)) %>%
    mutate(ID = paste(Chromatography, Polarity, RT, mz, sep = "_"),
           label = if_else(is.na(Identification) | Identification %in% c("NA",""),
                           ID, Identification)) %>%
    transmute(label, ID) %>%
    filter(ID %in% biomarqueur_dairy_fat)) %>%
  rename("Identification" = "label",
         "Biomarqueur" = "ID")

# Ajout des scores
identification <- rbind(identification,
                        data.frame(Identification = c("Milk multi-BFIs model", "Butter multi-BFIs model", "CFF dairy multi-BFIs model"),
                                   Biomarqueur = c("Score_lait", "Score_beurre", "Score_dairy_fat")))
# Mise en forme des données ----

# Création d'un vecteur avec tous les biomarqueurs
all_biomarkers <- c(biomarqueur_lait, biomarqueur_beurre, biomarqueur_dairy_fat)

# Selection des données au baseline des RCT + eMECA
data_ech <- subset(data_ech, Traitement == "Baseline" | Projet == "eMECA", select = c("Projet", "Sujet", "Projet_Sujet", all_biomarkers))

# Scaling
data_ech[,  all_biomarkers] <- scale(log(data_ech[,  all_biomarkers]))

# Merge avec les données de score multi-metabolites
data_ech <- merge(x = data_score, data_ech, by = c("Projet", "Sujet", "Projet_Sujet"))

# Merge avec les données de FFQ
data_ffq <- merge(x = data_ffq, y = data_ech, by = c("Projet", "Sujet"))

# Merge avec les données de R24W
data_r24w_q3 <- merge(x = data_r24w_q3, y = data_ech, by = c("Projet", "Sujet"))

data_r24w_mean <- merge(x = data_r24w_mean, y = data_ech, by = c("Projet", "Sujet"))

# Fonction pour régression linéaire multiple ----
reg_fun <- function(biomarqueur, data_df) {
  
  # Entraînement du modèle
  modele <- lm(as.formula(paste(biomarqueur, "~",
                                "Lait + Fromage + Yogourt + Creme + Beurre +
                                     Fruits + Legumes + Legumineuses + Grains_entiers + Grains_raffines + Viande + Poisson + Oeuf")),
               data = data_df)
  
  # Extraction des p-values
  results_p_value <- rownames_to_column(as.data.frame(summary(modele)$coefficients), var = "Aliment") %>%
    filter(Aliment != "(Intercept)") %>%
    rename("p_value" = "Pr(>|t|)") %>%
    select(Aliment, p_value) %>%
    mutate(p_sign = case_when(p_value < 0.0001 ~ "****",
                               p_value >= 0.0001 & p_value < 0.001 ~ "***",
                               p_value >= 0.001 & p_value < 0.01 ~ "**",
                               p_value >= 0.01 & p_value < 0.05 ~ "*",
                               p_value >= 0.05 ~ ""))
  
  # Extraction des r partiel
  results_corr_part <- rownames_to_column(as.data.frame(ols_correlations(modele)), var = "Aliment") %>%
    rename("r_partial" = "Partial") %>%
    select(Aliment, r_partial)
  
  # Graphique des postulats de régression
  print(ggarrange(ols_plot_resid_fit(modele), ols_plot_resid_qq(modele)))
  
  # Multicolinéarité (VIF)
  if(sum(vif(modele) > 10) > 0) {
    cat("!!! WARNING multicolinéarité (VIF > 10) !!!")
  }
  
  # Exportation des résultats
  return(full_join(x = results_corr_part, y = results_p_value, by = "Aliment") %>% mutate(Biomarqueur = biomarqueur))
}

# Analyse par régression linéaire ----

# Creation d'un vecteur avec les biomarqueurs et score a tester
var_reg <- c(biomarqueur_lait, "Score_lait", biomarqueur_beurre, "Score_beurre", biomarqueur_dairy_fat, "Score_dairy_fat")

# FFQ #

# Initialisation du data.frame
results_ffq <- data.frame(Aliment = NA, r_partial = NA, p_value = NA, p_sign = NA, Biomarqueur = NA)

# Boucle pour faire les regressions avec chacun des biomarqueurs
for(i in 1:length(var_reg)) {
  
  # Regression
  results <- reg_fun(var_reg[i], data_ffq)
  
  # Stockage des resultats
  results_ffq <- rbind(results_ffq, results)
  
  # Effacer le data.frame temporaire
  rm(results)
}

# Retrait de la ligne de NA
results_ffq <- results_ffq[-1, ]

# Mean of 24h recalls #

# Initialisation du data.frame
results_r24w_mean <- data.frame(Aliment = NA, r_partial = NA, p_value = NA, p_sign = NA, Biomarqueur = NA)

# Boucle pour faire les regressions avec chacun des biomarqueurs
for(i in 1:length(var_reg)) {
  
  # Regression
  results <- reg_fun(var_reg[i], data_r24w_mean)
  
  # Stockage des resultats
  results_r24w_mean <- rbind(results_r24w_mean, results)
  
  # Effacer le data.frame temporaire
  rm(results)
}

# Retrait de la ligne de NA
results_r24w_mean <- results_r24w_mean[-1, ]

# Day-of-draw 24h recall #

# Initialisation du data.frame
results_r24w_q3 <- data.frame(Aliment = NA, r_partial = NA, p_value = NA, p_sign = NA, Biomarqueur = NA)

# Boucle pour faire les regressions avec chacun des biomarqueurs
for(i in 1:length(var_reg)) {
  
  # Regression
  results <- reg_fun(var_reg[i], data_r24w_q3)
  
  # Stockage des resultats
  results_r24w_q3 <- rbind(results_r24w_q3, results)
  
  # Effacer le data.frame temporaire
  rm(results)
}

# Retrait de la ligne de NA
results_r24w_q3 <- results_r24w_q3[-1, ]

# Mise en forme des résultats ----

# FFQ #
results_ffq <- left_join(x = results_ffq %>% 
                           mutate(Aliment = case_match(Aliment,
                                                       "Lait" ~ "Milk",
                                                       "Fromage" ~ "Cheese",
                                                       "Yogourt" ~ "Yogurt",
                                                       "Creme" ~ "Cream",
                                                       "Beurre" ~ "Butter",
                                                       "Legumes" ~ "Vegetables",
                                                       "Legumineuses" ~ "Legumes",
                                                       "Grains_entiers" ~ "Whole grains",
                                                       "Grains_raffines" ~ "Refined grains",
                                                       "Viande" ~ "Meat",
                                                       "Poisson" ~ "Fish",
                                                       "Oeuf" ~ "Eggs",
                                                       .default = Aliment)),
                         y = identification,
                         by = "Biomarqueur")

# Mean of 24h recalls #
results_r24w_mean <- left_join(x = results_r24w_mean %>% 
                                 mutate(Aliment = case_match(Aliment,
                                                             "Lait" ~ "Milk",
                                                             "Fromage" ~ "Cheese",
                                                             "Yogourt" ~ "Yogurt",
                                                             "Creme" ~ "Cream",
                                                             "Beurre" ~ "Butter",
                                                             "Legumes" ~ "Vegetables",
                                                             "Legumineuses" ~ "Legumes",
                                                             "Grains_entiers" ~ "Whole grains",
                                                             "Grains_raffines" ~ "Refined grains",
                                                             "Viande" ~ "Meat",
                                                             "Poisson" ~ "Fish",
                                                             "Oeuf" ~ "Eggs",
                                                             .default = Aliment)),
                               y = identification,
                               by = "Biomarqueur")

# Day-of-draw 24h recall #
results_r24w_q3 <- left_join(x = results_r24w_q3 %>% 
                               mutate(Aliment = case_match(Aliment,
                                                           "Lait" ~ "Milk",
                                                           "Fromage" ~ "Cheese",
                                                           "Yogourt" ~ "Yogurt",
                                                           "Creme" ~ "Cream",
                                                           "Beurre" ~ "Butter",
                                                           "Legumes" ~ "Vegetables",
                                                           "Legumineuses" ~ "Legumes",
                                                           "Grains_entiers" ~ "Whole grains",
                                                           "Grains_raffines" ~ "Refined grains",
                                                           "Viande" ~ "Meat",
                                                           "Poisson" ~ "Fish",
                                                           "Oeuf" ~ "Eggs",
                                                           .default = Aliment)),
                             y = identification,
                             by = "Biomarqueur")

# FFQ #
ffq_sign <- results_ffq %>%
  select(Aliment, Identification, p_sign) %>%
  pivot_wider(names_from = "Identification", values_from = "p_sign") %>%
  column_to_rownames(var = "Aliment") %>% as.matrix()

graph_ffq <- Heatmap(results_ffq %>%
                       select(Aliment, Identification, r_partial) %>%
                       pivot_wider(names_from = "Identification", values_from = "r_partial") %>%
                       column_to_rownames(var = "Aliment") %>% as.matrix(),
                     name = "Partial ρ",
                     col = colorRamp2(c(-0.5, 0, 0.5), c("red", "white", "blue")),
                     rect_gp = gpar(col = "black"), # cell borders
                     show_row_dend = FALSE,
                     show_column_dend = FALSE,
                     cluster_rows = FALSE,
                     cluster_columns = FALSE,
                     row_names_side = "left",
                     column_names_side = "top",
                     column_names_gp = gpar(fontface = ifelse(colnames(results_ffq %>%
                                                                         select(Aliment, Identification, r_partial) %>%
                                                                         pivot_wider(names_from = "Identification", values_from = "r_partial") %>%
                                                                         column_to_rownames(var = "Aliment") %>%
                                                                         as.matrix()) %in% c("Milk multi-BFIs model", "Butter multi-BFIs model", "CFF dairy multi-BFIs model"), 
                                                              "bold", "plain")),
                     cell_fun = function(j, i, x, y, width, height, fill) {
                       grid.text(ffq_sign[i, j], x, y, gp = gpar(fontsize = 10.5, fontface = "bold"))
                     })

# Mean of 24h recalls #
r24w_mean_sign <- results_r24w_mean %>%
  select(Aliment, Identification, p_sign) %>%
  pivot_wider(names_from = "Identification", values_from = "p_sign") %>%
  column_to_rownames(var = "Aliment") %>% as.matrix()

graph_r24w_mean <- Heatmap(results_r24w_mean %>%
                       select(Aliment, Identification, r_partial) %>%
                       pivot_wider(names_from = "Identification", values_from = "r_partial") %>%
                       column_to_rownames(var = "Aliment") %>% as.matrix(),
                     name = "Partial ρ",
                     col = colorRamp2(c(-0.5, 0, 0.5), c("red", "white", "blue")),
                     rect_gp = gpar(col = "black"), # cell borders
                     show_row_dend = FALSE,
                     show_column_dend = FALSE,
                     cluster_rows = FALSE,
                     cluster_columns = FALSE,
                     row_names_side = "left",
                     column_names_side = "top",
                     column_names_gp = gpar(fontface = ifelse(colnames(results_r24w_mean %>%
                                                                         select(Aliment, Identification, r_partial) %>%
                                                                         pivot_wider(names_from = "Identification", values_from = "r_partial") %>%
                                                                         column_to_rownames(var = "Aliment") %>%
                                                                         as.matrix()) %in% c("Milk multi-BFIs model", "Butter multi-BFIs model", "CFF dairy multi-BFIs model",
                                                                                             "Any dairy multi-BFIs model"), 
                                                              "bold", "plain")),
                     cell_fun = function(j, i, x, y, width, height, fill) {
                       grid.text(r24w_mean_sign[i, j], x, y, gp = gpar(fontsize = 10.5, fontface = "bold"))
                     })


# Day-of-draw 24h recall #
r24w_q3_sign <- results_r24w_q3 %>%
  select(Aliment, Identification, p_sign) %>%
  pivot_wider(names_from = "Identification", values_from = "p_sign") %>%
  column_to_rownames(var = "Aliment") %>% as.matrix()

graph_r24w_q3 <- Heatmap(results_r24w_q3 %>%
                       select(Aliment, Identification, r_partial) %>%
                       pivot_wider(names_from = "Identification", values_from = "r_partial") %>%
                       column_to_rownames(var = "Aliment") %>% as.matrix(),
                     name = "Partial ρ",
                     col = colorRamp2(c(-0.5, 0, 0.5), c("red", "white", "blue")),
                     rect_gp = gpar(col = "black"), # cell borders
                     show_row_dend = FALSE,
                     show_column_dend = FALSE,
                     cluster_rows = FALSE,
                     cluster_columns = FALSE,
                     row_names_side = "left",
                     column_names_side = "top",
                     column_names_gp = gpar(fontface = ifelse(colnames(results_r24w_q3 %>%
                                                                         select(Aliment, Identification, r_partial) %>%
                                                                         pivot_wider(names_from = "Identification", values_from = "r_partial") %>%
                                                                         column_to_rownames(var = "Aliment") %>%
                                                                         as.matrix()) %in% c("Milk multi-BFIs model", "Butter multi-BFIs model", "CFF dairy multi-BFIs model",
                                                                                             "Any dairy multi-BFIs model"), 
                                                              "bold", "plain")),
                     cell_fun = function(j, i, x, y, width, height, fill) {
                       grid.text(r24w_q3_sign[i, j], x, y, gp = gpar(fontsize = 10.5, fontface = "bold"))
                     })

# Exportation des figures ----

# Figure 6
setwd("C:/Users/jackl/OneDrive/Bureau/École/Postdoctorat/Chercheur/Jean-Philippe/PLC/Article/Biomarqueurs/Figures")

tiff(filename = "Figure 6 - Heatmap correlation.tiff", width = 16, height = 19, units = "cm", res = 600, compression = "lzw")

draw(graph_ffq %v% graph_r24w_mean, ht_gap = unit(1, "cm"))

dev.off()

# Figure S4
tiff(filename = "Figure S4 - Heatmap correlation - Day-of-draw R24.tiff", width = 16, height = 13, units = "cm", res = 600, compression = "lzw")

draw(graph_r24w_q3)

dev.off()

