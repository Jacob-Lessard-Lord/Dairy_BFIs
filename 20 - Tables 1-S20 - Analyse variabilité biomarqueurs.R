# Date : 16 septembre 2025
# Auteur : Jacob Lessard-Lord

# Projet PLC

# Analyse des facteurs expliquant la variabilité des biomarqueurs

# Untargeted 

# Set working directory ----
setwd("C:/Users/jackl/OneDrive/Bureau/École/Postdoctorat/Chercheur/Jean-Philippe/PLC/Processing_PLC/Untargeted/Final_Processing/Data")

# Chargement des packages requis ----
library(tidyverse)
library(ggpubr)
library(car)
library(olsrr)
library(rsq)
library(readxl)
library(openxlsx)

# Importation des données métabolomiques ----
load("Metabo/Untargeted - ALL.rda")

# Importation  des biomarqueurs ----
load("Metabo/Biomarqueurs finaux.rda")

# Importation des données cliniques ----
load("FANI/Data clinique - RCT.rda")

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

# Mise en forme des données ----

# Création d'un vecteur avec tous les biomarqueurs
all_biomarkers <- c(biomarqueur_lait, biomarqueur_beurre, biomarqueur_dairy_fat)

# Selection des données post-dietes RCT
data_ech <- subset(data_ech, subset = Traitement %in% c("Lait", "Fromage", "Beurre", "Controle"),
                   select = c("Projet", "Sujet", "Projet_Sujet", "Traitement", all_biomarkers))

# Scaling
data_ech[,  all_biomarkers] <- scale(log(data_ech[,  all_biomarkers]))

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
                         data_ech, by = c("Projet", "Sujet", "Traitement")) %>%
  filter(!Projet_Sujet %in% c("HDL_24", "HDL_68"))

# Analyse de la variabilité ----
reg_fun <- function(biomarqueur, exposition, data_df) {
  
  # Entraînement du modèle
  modele <- lm(as.formula(paste(biomarqueur, "~ sex + age + bmi +", exposition)), data = data_df)
  
  # Ensure downstream functions (rsq.partial, drop1, update, etc.) can refit
  # without needing a 'data_df' object in the parent frame:
  modele$call$data <- data_df
  
  # Multicolinearite - les VIF doivent etre < 10 
  print(ols_vif_tol(modele))
  
  # Vérification des hypothèses du modèle
  print(ols_plot_resid_fit(modele)) # Graphique des résidus - Homogénéité
  print(ols_plot_resid_qq(modele)) # Q-Q Plot - Normalité des résidus
  print(ols_test_normality(modele)) # Test de normalité
  
  # Influence des observations
  print(influenceIndexPlot(modele))
  
  # Extraction des résultats
  results <- merge(x = rownames_to_column(as.data.frame(summary(modele)$coefficients), var = "Variable"),
                   y = rownames_to_column(as.data.frame(confint(modele)), var = "Variable"),
                   by = "Variable", sort = FALSE)[-1, ] # Retrait de l'intercept
  
  # Renommer les colonnes
  colnames(results) <- c("Variable", "Estimate", "Std_error", "t_value", "p_value", "lower_ci", "upper_ci")
  
  # Intégrer les r2 partiels aux résultats
  results <- cbind(results,
                   data.frame(Variable2 = rsq.partial(modele)$variable,
                              partial_r2 = rsq.partial(modele)$partial.rsq))
  
  # WARNINGS NORMAUX #
  
  # Mise en forme des résultats
  results <- results %>%
    mutate(Estimate = round(Estimate, digits = 3),
           lower_ci = round(lower_ci, digits = 3),
           upper_ci = round(upper_ci, digits = 3),
           beta_ci = paste0(Estimate, " (", lower_ci, ", ", upper_ci, ")"),
           p_value = round(p_value, digits = 4),
           partial_r2 = round((partial_r2 * 100), digits = 1))
  
  # Sélection des résultats pertinents
  results <- results %>% select(c("Variable", "beta_ci", "partial_r2", "p_value"))
  
  # Ajout du r2 du modèle
  results <- rbind(
    results,
    data.frame(Variable = "Model",
               beta_ci = NA,
               partial_r2 = round(ols_regress(modele)$rsq * 100, 1),
               p_value = round(ols_regress(modele)$p, 4))
  )
  
  # Arrondissement des p-values + Ajout des variables Exposition et Biomarqueur
  results <- results %>%
    mutate(p_value = case_when(p_value < 0.0001 ~ "<0.0001",
                               p_value >= 0.0001 & p_value < 0.001 ~ as.character(round(p_value, 4)),
                               p_value >= 0.001 & p_value < 0.01 ~ as.character(round(p_value, 3)),
                               p_value >= 0.01 ~ as.character(round(p_value, 2))),
           Exposition = exposition,
           Biomarqueur = biomarqueur) %>%
    select(Exposition, Biomarqueur, Variable, beta_ci, partial_r2, p_value)
  
  return(results)
  
}

# Lait
results_M1_lait <- reg_fun(biomarqueur_lait[1], "Milk_intake", data_rct_lait %>% mutate(Milk_intake = (Milk_intake / 100))) # Warning normal
results_M2_lait <- reg_fun(biomarqueur_lait[2], "Milk_intake", data_rct_lait %>% mutate(Milk_intake = (Milk_intake / 100))) # Warning normal

# Fromage
results_CFFD1_fromage <- reg_fun(biomarqueur_dairy_fat[1], "Cheese_intake", data_rct_fromage %>% mutate(Cheese_intake = (Cheese_intake / 10))) # Warning normal
results_CFFD2_fromage <- reg_fun(biomarqueur_dairy_fat[2], "Cheese_intake", data_rct_fromage %>% mutate(Cheese_intake = (Cheese_intake / 10))) # Warning normal
results_CFFD3_fromage <- reg_fun(biomarqueur_dairy_fat[3], "Cheese_intake", data_rct_fromage %>% mutate(Cheese_intake = (Cheese_intake / 10))) # Warning normal
results_CFFD4_fromage <- reg_fun(biomarqueur_dairy_fat[4], "Cheese_intake", data_rct_fromage %>% mutate(Cheese_intake = (Cheese_intake / 10))) # Warning normal
results_CFFD5_fromage <- reg_fun(biomarqueur_dairy_fat[5], "Cheese_intake", data_rct_fromage %>% mutate(Cheese_intake = (Cheese_intake / 10))) # Warning normal

# Beurre
results_B1_beurre <- reg_fun(biomarqueur_beurre[1], "Butter_intake", data_rct_beurre) # Warning normal
results_B2_beurre <- reg_fun(biomarqueur_beurre[2], "Butter_intake", data_rct_beurre) # Warning normal
results_CFFD1_beurre <- reg_fun(biomarqueur_dairy_fat[1], "Butter_intake", data_rct_beurre) # Warning normal
results_CFFD2_beurre <- reg_fun(biomarqueur_dairy_fat[2], "Butter_intake", data_rct_beurre) # Warning normal
results_CFFD3_beurre <- reg_fun(biomarqueur_dairy_fat[3], "Butter_intake", data_rct_beurre) # Warning normal
results_CFFD4_beurre <- reg_fun(biomarqueur_dairy_fat[4], "Butter_intake", data_rct_beurre) # Warning normal
results_CFFD5_beurre <- reg_fun(biomarqueur_dairy_fat[5], "Butter_intake", data_rct_beurre) # Warning normal


# Mise en forme finale et exportation des résultats ----
write.xlsx(left_join(x = rbind(results_M1_lait, results_M2_lait, 
                               results_CFFD1_fromage, results_CFFD2_fromage, results_CFFD3_fromage, results_CFFD4_fromage, results_CFFD5_fromage,
                               results_B1_beurre, results_B2_beurre,
                               results_CFFD1_beurre, results_CFFD2_beurre, results_CFFD3_beurre, results_CFFD4_beurre, results_CFFD5_beurre),
                     y = identification, by = "Biomarqueur") %>%
             mutate(Biomarqueur = Identification,
                    Variable = str_replace(Variable, "sexMale", "Sex (female vs male)"),
                    Variable = str_replace(Variable, "age", "Age (years)"),
                    Variable = str_replace(Variable, "bmi", "BMI (kg/m2)"),
                    Variable = str_replace(Variable, "_intake", " intake (g/day)"),
                    Exposition = str_remove(Exposition, "_intake")) %>%
             select(-Identification), 
           file = "Tables 1-S13 - Biomarker variability.xlsx")