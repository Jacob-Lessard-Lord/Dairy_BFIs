# Date : 13 août 2025
# Auteur : Jacob Lessard-Lord

# Projet PLC

# ÉValuation de la corrélation de la quantité de produits laitiers rapportés par FFQ et R24W

# Correlation de Pearson et Spearman

# Untargeted 

# Set working directory ----
setwd("C:/Users/jackl/OneDrive/Bureau/École/Postdoctorat/Chercheur/Jean-Philippe/PLC/Processing_PLC/Untargeted/Final_Processing/Data")

# Chargement des packages requis ----
library(tidyverse)
library(psych)
library(ComplexHeatmap)
library(circlize)

# Importation des données nutritionnelles ----
load("FFQ/Data FFQ.rda")
load("R24W/Data R24W.rda")

# Importation des données métabolomiques ----
load("Metabo/Untargeted - ALL.rda")

# Mise en forme des données ----

# Création d'un data.frame avec seulement les données FFQ d'eMECA
data_ffq <- data_ffq %>% filter(Projet == "eMECA")

# Création d'un vecteur contenant tous les sujets qui ont un FFQ valide + 3 R24W
sujet_unique <- intersect(data_ffq$Sujet, data_r24w_q3$Sujet)

# Selection de ces sujets dans les données
data_ffq <- data_ffq %>% filter(Sujet %in% sujet_unique)
data_r24w_mean <- data_r24w_mean %>% filter(Sujet %in% sujet_unique)
data_r24w_q3 <- data_r24w_q3 %>% filter(Sujet %in% sujet_unique)

# Selection des données au baseline des RCT + eMECA
data_ech <- subset(data_ech, Traitement == "Baseline" | Projet == "eMECA", select = c("Projet", "Sujet", "Projet_Sujet"))

# Merge avec les données de FFQ
data_ffq <- merge(x = data_ffq, y = data_ech, by = c("Projet", "Sujet"))

# Merge avec les données de R24W
data_r24w_q3 <- merge(x = data_r24w_q3, y = data_ech, by = c("Projet", "Sujet"))
data_r24w_mean <- merge(x = data_r24w_mean, y = data_ech, by = c("Projet", "Sujet"))

rm(data_ech, metabolite, sujet_unique)

# Recodage du nom des variables

## FFQ
data_ffq <- data_ffq %>%
  select(Sujet, Dairy_fat, Lait, Fromage, Yogourt, Beurre, Creme,
         Fruits, Legumes, Legumineuses, Grains_entiers, Grains_raffines, Viande, Poisson, Oeuf) %>%
  rename("CFF dairy" = Dairy_fat,
         "Milk" = Lait,
         "Cheese" = Fromage,
         "Yogurt" = Yogourt,
         "Butter" = Beurre,
         "Cream" = Creme,
         "Fruits" = Fruits,
         "Vegetables" = Legumes,
         "Legumes" = Legumineuses,
         "Whole grains" = Grains_entiers,
         "Refined grains" = Grains_raffines,
         "Meat" = Viande,
         "Fish" = Poisson,
         "Eggs" = Oeuf) %>%
  rename_with(~ paste0(.x, " - FFQ"), -Sujet)

## R24W - Q3
data_r24w_q3 <- data_r24w_q3 %>%
  select(Sujet, Dairy_fat, Lait, Fromage, Yogourt, Beurre, Creme,
         Fruits, Legumes, Legumineuses, Grains_entiers, Grains_raffines, Viande, Poisson, Oeuf) %>%
  rename("CFF dairy" = Dairy_fat,
         "Milk" = Lait,
         "Cheese" = Fromage,
         "Yogurt" = Yogourt,
         "Butter" = Beurre,
         "Cream" = Creme,
         "Fruits" = Fruits,
         "Vegetables" = Legumes,
         "Legumes" = Legumineuses,
         "Whole grains" = Grains_entiers,
         "Refined grains" = Grains_raffines,
         "Meat" = Viande,
         "Fish" = Poisson,
         "Eggs" = Oeuf) %>%
  rename_with(~ paste0(.x, " - R24"), -Sujet)

## Mean of 3x R24W
data_r24w_mean <- data_r24w_mean %>%
  select(Sujet, Dairy_fat, Lait, Fromage, Yogourt, Beurre, Creme,
         Fruits, Legumes, Legumineuses, Grains_entiers, Grains_raffines, Viande, Poisson, Oeuf) %>%
  rename("CFF dairy" = Dairy_fat,
         "Milk" = Lait,
         "Cheese" = Fromage,
         "Yogurt" = Yogourt,
         "Butter" = Beurre,
         "Cream" = Creme,
         "Fruits" = Fruits,
         "Vegetables" = Legumes,
         "Legumes" = Legumineuses,
         "Whole grains" = Grains_entiers,
         "Refined grains" = Grains_raffines,
         "Meat" = Viande,
         "Fish" = Poisson,
         "Eggs" = Oeuf) %>%
  rename_with(~ paste0(.x, " - mR24"), -Sujet)

# Merge des 3 datasets
data_self_report <- full_join(x = full_join(x = data_ffq %>% mutate(Sujet = factor(Sujet)),
                                            y = data_r24w_q3,
                                            by = "Sujet"),
                              y = data_r24w_mean,
                              by = "Sujet")

rm(data_ffq, data_r24w_q3, data_r24w_mean)

# Changer l'ordre des colonnes
data_self_report <- data_self_report %>%
  select("Sujet",
         "Milk - FFQ", "Milk - R24", "Milk - mR24",
         "Yogurt - FFQ", "Yogurt - R24", "Yogurt - mR24",
         "Cheese - FFQ", "Cheese - R24", "Cheese - mR24",
         "Cream - FFQ", "Cream - R24", "Cream - mR24",
         "Butter - FFQ", "Butter - R24", "Butter - mR24",
         "CFF dairy - FFQ", "CFF dairy - R24", "CFF dairy - mR24",
         "Fruits - FFQ", "Fruits - R24", "Fruits - mR24",
         "Vegetables - FFQ", "Vegetables - R24", "Vegetables - mR24",
         "Legumes - FFQ", "Legumes - R24", "Legumes - mR24",
         "Whole grains - FFQ", "Whole grains - R24", "Whole grains - mR24",
         "Refined grains - FFQ", "Refined grains - R24", "Refined grains - mR24",
         "Meat - FFQ", "Meat - R24", "Meat - mR24",
         "Fish - FFQ", "Fish - R24", "Fish - mR24",
         "Eggs - FFQ", "Eggs - R24", "Eggs - mR24")

# Graphique - Dairy only ----

# Correlation de Spearman
corr_self_report <- corr.test(data_self_report %>% select(c("Milk - FFQ", "Milk - R24", "Milk - mR24",
                                                            "Yogurt - FFQ", "Yogurt - R24", "Yogurt - mR24",
                                                            "Cheese - FFQ", "Cheese - R24", "Cheese - mR24",
                                                            "Cream - FFQ", "Cream - R24", "Cream - mR24",
                                                            "Butter - FFQ", "Butter - R24", "Butter - mR24",
                                                            "CFF dairy - FFQ", "CFF dairy - R24", "CFF dairy - mR24")),
                              data_self_report %>% select(c("Milk - FFQ", "Milk - R24", "Milk - mR24",
                                                            "Yogurt - FFQ", "Yogurt - R24", "Yogurt - mR24",
                                                            "Cheese - FFQ", "Cheese - R24", "Cheese - mR24",
                                                            "Cream - FFQ", "Cream - R24", "Cream - mR24",
                                                            "Butter - FFQ", "Butter - R24", "Butter - mR24",
                                                            "CFF dairy - FFQ", "CFF dairy - R24", "CFF dairy - mR24")),
                              method = "spearman")

# Matrice significativité
mat_sign_self_report <- matrix("", nrow = nrow(corr_self_report$p), ncol = ncol(corr_self_report$p),
                               dimnames = dimnames(corr_self_report$p))

mat_sign_self_report[corr_self_report$p < 0.0001] <- "****"
mat_sign_self_report[corr_self_report$p >= 0.0001 & corr_self_report$p < 0.001] <- "***"
mat_sign_self_report[corr_self_report$p >= 0.001 & corr_self_report$p < 0.01] <- "**"
mat_sign_self_report[corr_self_report$p >= 0.01 & corr_self_report$p < 0.05] <- "*"

# Exportation de la figure
setwd("C:/Users/jackl/OneDrive/Bureau/École/Postdoctorat/Chercheur/Jean-Philippe/PLC/Article/Biomarqueurs/Figures")

tiff(filename = "Figure S5 - Correlation - FFQ - R24W - Dairy.tiff", width = 20, height = 17.5, units = "cm", res = 600, compression = "lzw")

Heatmap(corr_self_report$r,
        name = "Spearman ρ",
        col = colorRamp2(c(-0.5, 0, 0.5, 1), c("red", "white", "blue", "darkblue")),
        na_col = "white",                 # cache la moitié supérieure
        rect_gp = gpar(col = "black"), # cell borders
        show_row_dend = FALSE,
        show_column_dend = FALSE,
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        row_names_side = "left",
        column_names_side = "top",
        #column_names_rot = 45,
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(mat_sign_self_report[i, j], x, y, gp = gpar(fontsize = 10.5, fontface = "bold"))
        })

dev.off()

# Graphique - Dairy vs Other ----


# Correlation de Spearman
corr_self_report <- corr.test(data_self_report %>% select(c("Fruits - FFQ", "Fruits - R24", "Fruits - mR24",
                                                            "Vegetables - FFQ", "Vegetables - R24", "Vegetables - mR24",
                                                            "Legumes - FFQ", "Legumes - R24", "Legumes - mR24",
                                                            "Whole grains - FFQ", "Whole grains - R24", "Whole grains - mR24",
                                                            "Refined grains - FFQ", "Refined grains - R24", "Refined grains - mR24",
                                                            "Meat - FFQ", "Meat - R24", "Meat - mR24",
                                                            "Fish - FFQ", "Fish - R24", "Fish - mR24",
                                                            "Eggs - FFQ", "Eggs - R24", "Eggs - mR24")),
                              data_self_report %>% select(c("Milk - FFQ", "Milk - R24", "Milk - mR24",
                                                            "Yogurt - FFQ", "Yogurt - R24", "Yogurt - mR24",
                                                            "Cheese - FFQ", "Cheese - R24", "Cheese - mR24",
                                                            "Cream - FFQ", "Cream - R24", "Cream - mR24",
                                                            "Butter - FFQ", "Butter - R24", "Butter - mR24",
                                                            "CFF dairy - FFQ", "CFF dairy - R24", "CFF dairy - mR24")),
                              method = "spearman")

# Matrice significativité
mat_sign_self_report <- matrix("", nrow = nrow(corr_self_report$p), ncol = ncol(corr_self_report$p),
                               dimnames = dimnames(corr_self_report$p))

mat_sign_self_report[corr_self_report$p < 0.0001] <- "****"
mat_sign_self_report[corr_self_report$p >= 0.0001 & corr_self_report$p < 0.001] <- "***"
mat_sign_self_report[corr_self_report$p >= 0.001 & corr_self_report$p < 0.01] <- "**"
mat_sign_self_report[corr_self_report$p >= 0.01 & corr_self_report$p < 0.05] <- "*"

# Exportation de la figure
setwd("C:/Users/jackl/OneDrive/Bureau/École/Postdoctorat/Chercheur/Jean-Philippe/PLC/Article/Biomarqueurs/Figures")

tiff(filename = "Figure S6 - Correlation - FFQ - R24W - Dairy vs Other.tiff", width = 20, height = 20, units = "cm", res = 600, compression = "lzw")

Heatmap(corr_self_report$r,
        name = "Spearman ρ",
        col = colorRamp2(c(-0.5, 0, 0.5, 1), c("red", "white", "blue", "darkblue")),
        na_col = "white",                 # cache la moitié supérieure
        rect_gp = gpar(col = "black"), # cell borders
        show_row_dend = FALSE,
        show_column_dend = FALSE,
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        row_names_side = "left",
        column_names_side = "top",
        #column_names_rot = 45,
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(mat_sign_self_report[i, j], x, y, gp = gpar(fontsize = 10.5, fontface = "bold"))
        })

dev.off()

