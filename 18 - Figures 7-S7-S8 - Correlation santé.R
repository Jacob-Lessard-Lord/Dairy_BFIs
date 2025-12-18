# Date : 15 septembre 2025
# Auteur : Jacob Lessard-Lord

# Projet PLC

# Corrélation entre les signatures métabolomiques et les facteurs de risque

# Untargeted 

# Set working directory ----
setwd("C:/Users/jackl/OneDrive/Bureau/École/Postdoctorat/Chercheur/Jean-Philippe/PLC/Processing_PLC/Untargeted/Final_Processing/Data")

# Chargement des packages requis ----
library(tidyverse)
library(readxl)
library(ggpubr)
library(psych)

# Importation des données cliniques ----
load("FANI/Data clinique - eMECA.rda")

# Importation des données nutritionnelles ----
load("FFQ/Data FFQ.rda")
load("R24W/Data R24W.rda")

# Importation des données métabolomiques ----
load("Metabo/Untargeted - ALL.rda")

# Importation  des biomarqueurs ----
load("Metabo/Biomarqueurs finaux.rda")

# Importation des scores multi-métabolites ----
load("Metabo/Score multi-metabolite.rda")

# Mise en forme des données ----

# Création d'un vecteur avec tous les biomarqueurs
all_biomarkers <- c(biomarqueur_lait, biomarqueur_beurre, biomarqueur_dairy_fat)

# Selection des données au baseline des RCT + eMECA
data_ech <- subset(data_ech, Projet == "eMECA", select = c("Projet", "Sujet", "Projet_Sujet", all_biomarkers))

# Scaling
data_ech[,  all_biomarkers] <- scale(log(data_ech[,  all_biomarkers]))

# Merge avec les données de score multi-metabolites
data_ech <- merge(x = data_score, data_ech, by = c("Projet", "Sujet", "Projet_Sujet"))

# Merge avec les données cliniques
data_clinique <- merge(x = data_clinique, y = data_ech, by = c("Projet", "Sujet")) %>%
  mutate(sex = ifelse(test = sex == "Female", yes = 0, no = 1))

# Merge des données de FFQ et de R24w
data_nutrition <- merge(x = data_ffq %>% select(c("Projet", "Sujet", "Lait", "Beurre", "Dairy_fat", "Total_dairy")) %>%
                          rename_with(~ paste0(.x, "_FFQ"), -c(Projet, Sujet)),
                        y = data_r24w_mean %>% select(c("Projet", "Sujet", "Lait", "Beurre", "Dairy_fat", "Total_dairy")) %>%
                          rename_with(~ paste0(.x, "_mR24"), -c(Projet, Sujet)),
                        by = c("Projet", "Sujet"))

# Merge des données nutritionnelles avec les autres
data_clinique <- merge(x = data_clinique,
                       y = data_nutrition,
                       by = c("Projet", "Sujet"))

# Recodage des noms  des métabolites ----

# Importation des identifications

# Lait #
id_lait <- read_excel("Metabo/ID Biomarqueurs.xlsx", sheet = "Lait") %>% select(c(Chromatography, Polarity, RT, mz, Identification)) %>%
  mutate(ID = paste(Chromatography, Polarity, RT, mz, sep = "_"),
         label = if_else(is.na(Identification) | Identification %in% c("NA",""),
                         ID, Identification)) %>%
  transmute(label, ID) %>%
  deframe()

# Beurre #
id_beurre <- read_excel("Metabo/ID Biomarqueurs.xlsx", sheet = "Beurre") %>% select(c(Chromatography, Polarity, RT, mz, Identification)) %>%
  mutate(ID = paste(Chromatography, Polarity, RT, mz, sep = "_"),
         label = if_else(is.na(Identification) | Identification %in% c("NA",""),
                         ID, Identification)) %>%
  transmute(label, ID) %>%
  deframe()

# Dairy fat #
id_dairy_fat <- read_excel("Metabo/ID Biomarqueurs.xlsx", sheet = "Dairy_fat") %>% select(c(Chromatography, Polarity, RT, mz, Identification)) %>%
  mutate(ID = paste(Chromatography, Polarity, RT, mz, sep = "_"),
         label = if_else(is.na(Identification) | Identification %in% c("NA",""),
                         ID, Identification)) %>%
  transmute(label, ID) %>%
  deframe()

# Recodage des noms

# Lait #
biomarqueur_lait <- id_lait[id_lait %in% biomarqueur_lait]

# Beurre #
biomarqueur_beurre <- id_beurre[id_beurre %in% biomarqueur_beurre]

# Dairy fat #
biomarqueur_dairy_fat <- id_dairy_fat[id_dairy_fat %in% biomarqueur_dairy_fat]

# Création des vecteurs contenant les variables à corréler ----

# Variables d'exposition
variables_lait <- c("Lait_FFQ", "Lait_mR24", "Score_lait", biomarqueur_lait)
variables_beurre <- c("Beurre_FFQ", "Beurre_mR24", "Score_beurre", biomarqueur_beurre)
variables_dairy_fat <- c("Dairy_fat_FFQ", "Dairy_fat_mR24", "Score_dairy_fat", biomarqueur_dairy_fat)

# Variables santé
variables_sante <- c("waistc", "bmi", "chol", "ldlc", "hdlc", "nhdlc", "tg", "crp_neph", "glufast", "insfast2018", "hba1c", "sbp", "dbp")

# Covariables
covariables <- c("age", "sex")

# Séparation du jeu de données par exposition ---

# Lait
data_lait <- data_clinique %>%
  select(all_of(c(unname(variables_lait), variables_sante, covariables)))

# Beurre
data_beurre <- data_clinique %>%
  select(all_of(c(unname(variables_beurre), variables_sante, covariables)))

# Dairy fat
data_dairy_fat <- data_clinique %>%
  select(all_of(c(unname(variables_dairy_fat), variables_sante, covariables)))

# Corrélation et mise en forme ----

# Lait # 

# 1) Faire la corrélation de Spearman sur toutes les variables
correlation_lait_r <- cor(data_lait, method = "spearman")

# 2) Obtenir les corrélations partielles en ajustant pour l'âge et le sexe
correlation_lait_r <- partial.r(correlation_lait_r,
                                x = match(c(variables_lait, variables_sante), colnames(correlation_lait_r)),
                                y = match(covariables, colnames(correlation_lait_r)))

# 3) Selectionner seulement les corrélations entre l'exposition et les facteurs de risque
correlation_lait_r <- correlation_lait_r[variables_lait, variables_sante]

# 4) Obtenir les p-values
correlation_lait_p <- corr.p(correlation_lait_r, n = (nrow(data_lait) - 2), adjust = "none")

# 5) Mise en forme en data.frame
correlation_lait_r <- data.frame(correlation_lait_r) %>%
  rownames_to_column(var = "Nutrition") %>%
  pivot_longer(cols = -Nutrition, names_to = "Risk_factor", values_to = "Spearman_r")

correlation_lait_p <- data.frame(correlation_lait_p$p) %>%
  rownames_to_column(var = "Nutrition") %>%
  pivot_longer(cols = -Nutrition, names_to = "Risk_factor", values_to = "p_value")

# 6) Merge des 2 data.frames
correlation_lait <- full_join(x = correlation_lait_r, y = correlation_lait_p, by = c("Nutrition", "Risk_factor"))

rm(correlation_lait_r, correlation_lait_p)

# 7) Mise en forme finale du nom des variables
correlation_lait <- correlation_lait %>%
  mutate(`Assessment method` = factor(case_when(str_detect(Nutrition, "FFQ") ~ "FFQ",
                                                str_detect(Nutrition, "mR24") ~ "mR24",
                                                str_detect(Nutrition, "Score") ~ "Multi-BFIs model",
                                                str_detect(Nutrition, "RPLC_Pos_0.75_160.1302mz") ~ "δ-Valerobetaine",
                                                str_detect(Nutrition, "RPLC_Pos_1.07_158.1170mz") ~ "Homostachydrine"),
                                      levels = c("FFQ", "mR24", "δ-Valerobetaine", "Homostachydrine", "Multi-BFIs model")),
         Risk_factor = factor(case_match(Risk_factor,
                                         "waistc" ~ "Waist circumference",
                                         "bmi" ~ "BMI",
                                         "chol" ~ "Total cholesterol",
                                         "ldlc" ~ "LDL-C",
                                         "hdlc" ~ "HDL-C",
                                         "nhdlc" ~ "Non-HDL-C",
                                         "tg" ~ "Triglycerides",
                                         "crp_neph" ~ "CRP",
                                         "glufast" ~ "Glucose",
                                         "insfast2018" ~ "Insulin",
                                         "hba1c" ~ "HbA1c",
                                         "sbp" ~ "Systolic BP",
                                         "dbp" ~ "Diastolic BP"),
                              levels = c("Diastolic BP", "Systolic BP", "HbA1c", "Glucose", "Insulin", "CRP",
                                         "Triglycerides", "Non-HDL-C", "HDL-C", "LDL-C", "Total cholesterol", "BMI","Waist circumference")),
         `Statistical\nsignificativity` = factor(case_when(p_value < 0.05 ~ "P-value < 0.05",
                                                           p_value >= 0.05 & p_value < 0.1 ~ "0.05 ≤ P-value < 0.10",
                                                           p_value > 0.1 ~ "P-value ≥ 0.10"),
                                                 levels = c("P-value < 0.05", "0.05 ≤ P-value < 0.10", "P-value ≥ 0.10")))



# Beurre # 

# 1) Faire la corrélation de Spearman sur toutes les variables
correlation_beurre_r <- cor(data_beurre, method = "spearman")

# 2) Obtenir les corrélations partielles en ajustant pour l'âge et le sexe
correlation_beurre_r <- partial.r(correlation_beurre_r,
                                  x = match(c(variables_beurre, variables_sante), colnames(correlation_beurre_r)),
                                  y = match(covariables, colnames(correlation_beurre_r)))

# 3) Selectionner seulement les corrélations entre l'exposition et les facteurs de risque
correlation_beurre_r <- correlation_beurre_r[variables_beurre, variables_sante]

# 4) Obtenir les p-values
correlation_beurre_p <- corr.p(correlation_beurre_r, n = (nrow(data_beurre) - 2), adjust = "none")

# 5) Mise en forme en data.frame
correlation_beurre_r <- data.frame(correlation_beurre_r) %>%
  rownames_to_column(var = "Nutrition") %>%
  pivot_longer(cols = -Nutrition, names_to = "Risk_factor", values_to = "Spearman_r")

correlation_beurre_p <- data.frame(correlation_beurre_p$p) %>%
  rownames_to_column(var = "Nutrition") %>%
  pivot_longer(cols = -Nutrition, names_to = "Risk_factor", values_to = "p_value")

# 6) Merge des 2 data.frames
correlation_beurre <- full_join(x = correlation_beurre_r, y = correlation_beurre_p, by = c("Nutrition", "Risk_factor"))

rm(correlation_beurre_r, correlation_beurre_p)

# 7) Mise en forme finale du nom des variables
correlation_beurre <- correlation_beurre %>%
  mutate(`Assessment method` = factor(case_when(str_detect(Nutrition, "FFQ") ~ "FFQ",
                                                str_detect(Nutrition, "mR24") ~ "mR24",
                                                str_detect(Nutrition, "Score") ~ "Multi-BFIs model",
                                                str_detect(Nutrition, "Lipido_Neg_14.665_736.52845mz") ~ "PE (P-17:0/20:4, O-17:1/20:4)",
                                                str_detect(Nutrition, "Lipido_Neg_14.722_824.5808mz") ~ "PC (P-17:0/20:4, O-17:1/20:4)"),
                                      levels = c("FFQ", "mR24", "PE (P-17:0/20:4, O-17:1/20:4)", "PC (P-17:0/20:4, O-17:1/20:4)", "Multi-BFIs model")),
         Risk_factor = factor(case_match(Risk_factor,
                                         "waistc" ~ "Waist circumference",
                                         "bmi" ~ "BMI",
                                         "chol" ~ "Total cholesterol",
                                         "ldlc" ~ "LDL-C",
                                         "hdlc" ~ "HDL-C",
                                         "nhdlc" ~ "Non-HDL-C",
                                         "tg" ~ "Triglycerides",
                                         "crp_neph" ~ "CRP",
                                         "glufast" ~ "Glucose",
                                         "insfast2018" ~ "Insulin",
                                         "hba1c" ~ "HbA1c",
                                         "sbp" ~ "Systolic BP",
                                         "dbp" ~ "Diastolic BP"),
                              levels = c("Diastolic BP", "Systolic BP", "HbA1c", "Glucose", "Insulin", "CRP",
                                         "Triglycerides", "Non-HDL-C", "HDL-C", "LDL-C", "Total cholesterol", "BMI","Waist circumference")),
         `Statistical\nsignificativity` = factor(case_when(p_value < 0.05 ~ "P-value < 0.05",
                                                           p_value >= 0.05 & p_value < 0.1 ~ "0.05 ≤ P-value < 0.10",
                                                           p_value > 0.1 ~ "P-value ≥ 0.10"),
                                                 levels = c("P-value < 0.05", "0.05 ≤ P-value < 0.10", "P-value ≥ 0.10")))


# Dairy fat # 

# 1) Faire la corrélation de Spearman sur toutes les variables
correlation_dairy_fat_r <- cor(data_dairy_fat, method = "spearman")

# 2) Obtenir les corrélations partielles en ajustant pour l'âge et le sexe
correlation_dairy_fat_r <- partial.r(correlation_dairy_fat_r,
                                     x = match(c(variables_dairy_fat, variables_sante), colnames(correlation_dairy_fat_r)),
                                     y = match(covariables, colnames(correlation_dairy_fat_r)))

# 3) Selectionner seulement les corrélations entre l'exposition et les facteurs de risque
correlation_dairy_fat_r <- correlation_dairy_fat_r[variables_dairy_fat, variables_sante]

# 4) Obtenir les p-values
correlation_dairy_fat_p <- corr.p(correlation_dairy_fat_r, n = (nrow(data_dairy_fat) - 2), adjust = "none")

# 5) Mise en forme en data.frame
correlation_dairy_fat_r <- data.frame(correlation_dairy_fat_r) %>%
  rownames_to_column(var = "Nutrition") %>%
  pivot_longer(cols = -Nutrition, names_to = "Risk_factor", values_to = "Spearman_r")

correlation_dairy_fat_p <- data.frame(correlation_dairy_fat_p$p) %>%
  rownames_to_column(var = "Nutrition") %>%
  pivot_longer(cols = -Nutrition, names_to = "Risk_factor", values_to = "p_value")

# 6) Merge des 2 data.frames
correlation_dairy_fat <- full_join(x = correlation_dairy_fat_r, y = correlation_dairy_fat_p, by = c("Nutrition", "Risk_factor"))

rm(correlation_dairy_fat_r, correlation_dairy_fat_p)

# 7) Mise en forme finale du nom des variables
correlation_dairy_fat <- correlation_dairy_fat %>%
  mutate(`Assessment method` = factor(case_when(str_detect(Nutrition, "FFQ") ~ "FFQ",
                                                str_detect(Nutrition, "mR24") ~ "mR24",
                                                str_detect(Nutrition, "Score") ~ "Multi-BFIs model",
                                                str_detect(Nutrition, "Lipido_Neg_13.982_761.58122mz") ~ "SM d18:1/17:0",
                                                str_detect(Nutrition, "Lipido_Neg_14.724_818.59093mz") ~ "PC 17:0/18:1",
                                                str_detect(Nutrition, "Lipido_Neg_15.626_871.68958mz") ~ "SM 43:2",
                                                str_detect(Nutrition, "Lipido_Neg_16.06_873.70541mz") ~ "SM 43:1",
                                                str_detect(Nutrition, "Lipido_Neg_16.155_698.62159mz") ~ "Cer d18:1/25:0"),
                                      levels = c("FFQ", "mR24", "SM d18:1/17:0", "SM 43:1", "Cer d18:1/25:0", "SM 43:2", "PC 17:0/18:1", "Multi-BFIs model")),
         Risk_factor = factor(case_match(Risk_factor,
                                         "waistc" ~ "Waist circumference",
                                         "bmi" ~ "BMI",
                                         "chol" ~ "Total cholesterol",
                                         "ldlc" ~ "LDL-C",
                                         "hdlc" ~ "HDL-C",
                                         "nhdlc" ~ "Non-HDL-C",
                                         "tg" ~ "Triglycerides",
                                         "crp_neph" ~ "CRP",
                                         "glufast" ~ "Glucose",
                                         "insfast2018" ~ "Insulin",
                                         "hba1c" ~ "HbA1c",
                                         "sbp" ~ "Systolic BP",
                                         "dbp" ~ "Diastolic BP"),
                              levels = c("Diastolic BP", "Systolic BP", "HbA1c", "Glucose", "Insulin", "CRP",
                                         "Triglycerides", "Non-HDL-C", "HDL-C", "LDL-C", "Total cholesterol", "BMI","Waist circumference")),
         `Statistical\nsignificativity` = factor(case_when(p_value < 0.05 ~ "P-value < 0.05",
                                                           p_value >= 0.05 & p_value < 0.1 ~ "0.05 ≤ P-value < 0.10",
                                                           p_value > 0.1 ~ "P-value ≥ 0.10"),
                                                 levels = c("P-value < 0.05", "0.05 ≤ P-value < 0.10", "P-value ≥ 0.10")))

# Graphique ----
setwd("C:/Users/jackl/OneDrive/Bureau/École/Postdoctorat/Chercheur/Jean-Philippe/PLC/Article/Biomarqueurs/Figures")

# Lait
tiff(filename = "Figure 7 - Correlation risk factor - Milk.tiff", width = 20, height = 25, units = "cm", res = 600, compression = "lzw")

ggplot(correlation_lait,
       aes(x = Risk_factor, y = Spearman_r,
           colour = `Assessment method`,
           shape  = `Statistical\nsignificativity`,
           group  = `Assessment method`)) +
  geom_linerange(aes(ymin = pmin(0, Spearman_r),
                     ymax = pmax(0, Spearman_r)),
                 position = position_dodge2(width = .75, preserve = "single"), linetype = "dotted", linewidth = 1) +
  geom_point(position = position_dodge2(width = .75, preserve = "single"), size = 2.2) +
  scale_colour_manual(values = c("#0077BB", "#EE7733", "#228833", "#CC3311", "#009988")) +
  scale_shape_manual(values = c(19, 10, 1)) +
  geom_hline(yintercept = 0) +
  geom_hline(yintercept = c(-0.1, 0.1, 0.2, 0.3), colour = "grey50", linetype = "dashed") +
  coord_flip() +
  labs(title = "Milk", y = "Spearman correlation coefficient") +
  theme_pubr() +
  guides(shape  = guide_legend(nrow = 1),
         colour = guide_legend(nrow = 2, byrow = TRUE)) +
  theme(axis.title.y = element_blank(),
        plot.title = element_text(hjust = .5, face = "bold"),
        legend.text = element_text(size = 11),
        legend.position = "bottom",
        legend.box      = "vertical",
        legend.spacing.y = unit(2, "mm"))

dev.off()


# Beurre
tiff(filename = "Figure S7 - Correlation risk factor - Butter.tiff", width = 21, height = 28, units = "cm", res = 600, compression = "lzw")

ggplot(correlation_beurre,
       aes(x = Risk_factor, y = Spearman_r,
           colour = `Assessment method`,
           shape  = `Statistical\nsignificativity`,
           group  = `Assessment method`)) +
  geom_linerange(aes(ymin = pmin(0, Spearman_r),
                     ymax = pmax(0, Spearman_r)),
                 position = position_dodge2(width = .75, preserve = "single"), linetype = "dotted", linewidth = 1) +
  geom_point(position = position_dodge2(width = .75, preserve = "single"), size = 2.2) +
  scale_colour_manual(values = c("#0077BB", "#EE7733", "#228833", "#CC3311", "#009988", "#EE3377")) +
  scale_shape_manual(values = c(19, 10, 1)) +
  geom_hline(yintercept = 0) +
  geom_hline(yintercept = c(-0.1, 0.1, 0.2, 0.3), colour = "grey50", linetype = "dashed") +
  coord_flip() +
  labs(title = "Butter", y = "Spearman correlation coefficient") +
  theme_pubr() +
  guides(shape  = guide_legend(nrow = 1),
         colour = guide_legend(nrow = 3, byrow = TRUE)) +
  theme(axis.title.y = element_blank(),
        plot.title = element_text(hjust = .5, face = "bold"),
        legend.text = element_text(size = 11),
        legend.position = "bottom",
        legend.box      = "vertical",
        legend.spacing.y = unit(2, "mm"))

dev.off()

# Dairy fat
tiff(filename = "Figure S8 - Correlation risk factor - Dairy fat.tiff", width = 25, height = 35, units = "cm", res = 600, compression = "lzw")

ggplot(correlation_dairy_fat,
       aes(x = Risk_factor, y = Spearman_r,
           colour = `Assessment method`,
           shape  = `Statistical\nsignificativity`,
           group  = `Assessment method`)) +
  geom_linerange(aes(ymin = pmin(0, Spearman_r),
                     ymax = pmax(0, Spearman_r)),
                 position = position_dodge2(width = .75, preserve = "single"), linetype = "dotted", linewidth = 1) +
  geom_point(position = position_dodge2(width = .75, preserve = "single"), size = 2.2) +
  scale_colour_manual(values = c("#0077BB", "#EE7733", "#228833", "#CC3311", "#009988", "#EE3377", "#F0E442", "#6F4C9B")) +
  scale_shape_manual(values = c(19, 10, 1)) +
  geom_hline(yintercept = 0) +
  geom_hline(yintercept = c(-0.2, -0.1, 0.1, 0.2, 0.3, 0.4), colour = "grey50", linetype = "dashed") +
  scale_y_continuous(breaks = c(-0.2, -0.1, 0, 0.1, 0.2, 0.3, 0.4)) +
  coord_flip() +
  labs(title = "Concentrated full-fat dairy", y = "Spearman correlation coefficient") +
  theme_pubr() +
  theme(axis.title.y = element_blank(),
        plot.title = element_text(hjust = .5, face = "bold"),
        legend.text = element_text(size = 11),
        legend.position = "bottom",
        legend.box = "vertical",
        legend.spacing.y = unit(2, "mm"))

dev.off()
