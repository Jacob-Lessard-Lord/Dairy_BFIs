# Date : 4 août 2025
# Auteur : Jacob Lessard-Lord

# Projet PLC

# Identification de la signature metabolomique du lait et du fromage 

# Wilcoxon test

# Untargeted

# Set working directory ----
setwd("C:/Users/jackl/OneDrive/Bureau/École/Postdoctorat/Chercheur/Jean-Philippe/PLC/Processing_PLC/Untargeted/Final_Processing/Data")

# Chargement des packages requis ----
library(tidyverse)
library(ggpubr)
library(rstatix)
library(openxlsx)

# Importation des données ----
load("Metabo/Untargeted - ALL.rda")

# Importation et mise en forme des biomarqueurs ----
load("Metabo/XGB - Selection variables.rda")

rm(results_xgb_select_lait, results_xgb_select_fromage, results_xgb_select_beurre)

# Subset des données pour le lait, le fromage et le beurre ----

# Merge des biomarqueurs potentiels pour le lait, le fromage et le beurre
metabolite_lait <- unique(c(metabolite_lait, metabolite_lait_fromage))
metabolite_fromage <- unique(c(metabolite_fromage, metabolite_lait_fromage, metabolite_fromage_beurre))
metabolite_beurre <- unique(c(metabolite_beurre, metabolite_fromage_beurre))

rm(metabolite_lait_fromage, metabolite_fromage_beurre)

# Création d'un vecteur avec les biomarqueurs potentiels
metabolite_dairy <- unique(c(metabolite_lait, metabolite_fromage, metabolite_beurre))

# Sélection des biomarqueurs potentiels seulement
data_ech <- data_ech[, c("Projet", "Sujet", "Projet_Sujet","Traitement", metabolite_dairy)]

# Normalisation
data_ech[,  metabolite_dairy] <- scale(log(data_ech[,  metabolite_dairy]))

# Lait #
data_lait <- subset(data_ech, subset = Projet %in% c("LAIT", "GABA2") & Traitement %in% c("Controle", "Lait") & Projet_Sujet != "LAIT_16") %>%
  mutate(Traitement = droplevels(Traitement))

# Fromage #
data_fromage <- subset(data_ech, subset = Projet %in% c("HDL", "GABA2") & Traitement %in% c("Controle", "Fromage")) %>%
  mutate(Traitement = droplevels(Traitement))

# Beurre #
data_beurre <- subset(data_ech, subset = Projet == "HDL" & Traitement %in% c("Controle", "Beurre")& !Projet_Sujet %in% c("HDL_24", "HDL_68")) %>%
  mutate(Traitement = droplevels(Traitement))

# Calcul du delta pour la sélection des métabolites ----

# Lait #

## Initialisation du data.frame pour les delta
data_lait_delta <- pivot_wider(data = data_lait[, c("Projet", "Sujet", "Traitement", metabolite_dairy[1])], names_from = "Traitement",
                               values_from = metabolite_dairy[1])[, c(1:2)]

# Boucle pour calcul du delta
for(i in 1:length(metabolite_dairy)) {
  df_delta <- pivot_wider(data = data_lait[, c("Projet", "Sujet", "Traitement", metabolite_dairy[i])], names_from = "Traitement",
                          values_from = metabolite_dairy[i])
  df_delta[, metabolite_dairy[i]] <- df_delta$Lait - df_delta$Controle
  data_lait_delta <- merge(x = data_lait_delta, y = df_delta[, c("Projet", "Sujet", metabolite_dairy[i])], by = c("Projet", "Sujet"))
  rm(df_delta)
}

## Fromage ##

## Initialisation du data.frame pour les delta
data_fromage_delta <- pivot_wider(data = data_fromage[, c("Projet", "Sujet", "Traitement", metabolite_dairy[1])], names_from = "Traitement",
                                  values_from = metabolite_dairy[1])[, c(1:2)]

# Boucle pour calcul du delta
for(i in 1:length(metabolite_dairy)) {
  df_delta <- pivot_wider(data = data_fromage[, c("Projet", "Sujet", "Traitement", metabolite_dairy[i])], names_from = "Traitement",
                          values_from = metabolite_dairy[i])
  df_delta[, metabolite_dairy[i]] <- df_delta$Fromage - df_delta$Controle
  data_fromage_delta <- merge(x = data_fromage_delta, y = df_delta[, c("Projet", "Sujet", metabolite_dairy[i])], by = c("Projet", "Sujet"))
  rm(df_delta)
}

## Beurre ##

## Initialisation du data.frame pour les delta
data_beurre_delta <- pivot_wider(data = data_beurre[, c("Projet", "Sujet", "Traitement", metabolite_dairy[1])], names_from = "Traitement",
                                 values_from = metabolite_dairy[1])[, c(1:2)]

# Boucle pour calcul du delta
for(i in 1:length(metabolite_dairy)) {
  df_delta <- pivot_wider(data = data_beurre[, c("Projet", "Sujet", "Traitement", metabolite_dairy[i])], names_from = "Traitement",
                          values_from = metabolite_dairy[i])
  df_delta[, metabolite_dairy[i]] <- df_delta$Beurre - df_delta$Controle
  data_beurre_delta <- merge(x = data_beurre_delta, y = df_delta[, c("Projet", "Sujet", metabolite_dairy[i])], by = c("Projet", "Sujet"))
  rm(df_delta)
}


# Mise en forme longue des données ----
data_lait_long <- pivot_longer(data_lait, cols = all_of(metabolite_dairy), names_to = "Metabolite", values_to = "Intensity")
data_fromage_long <- pivot_longer(data_fromage, cols = all_of(metabolite_dairy), names_to = "Metabolite", values_to = "Intensity")
data_beurre_long <- pivot_longer(data_beurre, cols = all_of(metabolite_dairy), names_to = "Metabolite", values_to = "Intensity")

# Mise en forme data PLI ----
data_pli <- subset(data_ech, subset = Projet == "PLI" & Traitement %in% c("Controle", "Produits laitiers"))
data_pli_long <- pivot_longer(data_pli, cols = all_of(metabolite_dairy), names_to = "Metabolite", values_to = "Intensity")

# Delta PLI 

## Initialisation du data.frame pour les delta
data_pli_delta <- pivot_wider(data = data_pli[, c("Projet", "Sujet", "Traitement", metabolite_dairy[1])], names_from = "Traitement",
                                 values_from = metabolite_dairy[1])[, c(1:2)]

# Boucle pour calcul du delta
for(i in 1:length(metabolite_dairy)) {
  df_delta <- pivot_wider(data = data_pli[, c("Projet", "Sujet", "Traitement", metabolite_dairy[i])], names_from = "Traitement",
                          values_from = metabolite_dairy[i])
  df_delta[, metabolite_dairy[i]] <- df_delta$`Produits laitiers` - df_delta$Controle
  data_pli_delta <- merge(x = data_pli_delta, y = df_delta[, c("Projet", "Sujet", metabolite_dairy[i])], by = c("Projet", "Sujet"))
  rm(df_delta)
}

# Initialisation des tableaux supplémentaires ----

# Lait
tableau_lait <- tibble(Candidate = metabolite_lait)

# Fromage
tableau_fromage <- tibble(Candidate = metabolite_fromage)

# Beurre
tableau_beurre <- tibble(Candidate = metabolite_beurre)

# Step 2 - Test statistique de la plausibilité ----

# Lait
test_sign_lait <- subset(data_lait_long, subset = Metabolite %in% metabolite_lait) %>%
  group_by(Metabolite) %>%
  wilcox_test(Intensity ~ Traitement, paired = TRUE, alternative = "less") %>%
  adjust_pvalue(method = "none") %>%
  add_significance()

# Fromage
test_sign_fromage <- subset(data_fromage_long, subset = Metabolite %in% metabolite_fromage) %>%
  group_by(Metabolite) %>%
  wilcox_test(Intensity ~ Traitement, paired = TRUE, alternative = "less") %>%
  adjust_pvalue(method = "none") %>%
  add_significance()

# Beurre
test_sign_beurre <- subset(data_beurre_long, subset = Metabolite %in% metabolite_beurre) %>%
  group_by(Metabolite) %>%
  wilcox_test(Intensity ~ Traitement, paired = TRUE, alternative = "less") %>%
  adjust_pvalue(method = "none") %>%
  add_significance()



# Extraction des résultats pour tableau supplémentaire

## Lait

### Moyenne (95% CI)
tableau_lait <- left_join(x = tableau_lait,
                          y = data_lait_delta %>%
                            select(all_of(metabolite_lait)) %>%
                            get_summary_stats(type = "mean_ci") %>%
                            mutate(mean = round(mean, 2),
                                   ci = round(ci, 2),
                                   mean_ci = paste0(mean, " (", (mean - ci), ", ", (mean + ci), ")")) %>%
                            select(variable, mean_ci) %>%
                            rename(Candidate = variable,
                                   Step2_mean_ci = mean_ci),
                          by = "Candidate")

### p-value
tableau_lait <- left_join(x = tableau_lait,
                          y = test_sign_lait %>%
                            select(Metabolite, p.adj) %>%
                            mutate(p.adj = case_when(
                              p.adj < 0.0001 ~ "<0.0001",
                              p.adj >= 0.0001 & p.adj < 0.001 ~ as.character(round(p.adj, 4)),
                              p.adj >= 0.001 & p.adj < 0.01  ~ as.character(round(p.adj, 3)),
                              p.adj >= 0.01              ~ as.character(round(p.adj, 2)),
                              .default = NA)) %>% 
                            rename(Candidate = Metabolite,
                                   Step2_p_adj = p.adj),
                          by = "Candidate")

## Fromage

### Moyenne (95% CI)
tableau_fromage <- left_join(x = tableau_fromage,
                             y = data_fromage_delta %>%
                               select(all_of(metabolite_fromage)) %>%
                               get_summary_stats(type = "mean_ci") %>%
                               mutate(mean = round(mean, 2),
                                      ci = round(ci, 2),
                                      mean_ci = paste0(mean, " (", (mean - ci), ", ", (mean + ci), ")")) %>%
                               select(variable, mean_ci) %>%
                               rename(Candidate = variable,
                                      Step2_mean_ci = mean_ci),
                             by = "Candidate")

### p-value
tableau_fromage <- left_join(x = tableau_fromage,
                             y = test_sign_fromage %>%
                               select(Metabolite, p.adj) %>%
                               mutate(p.adj = case_when(
                                 p.adj < 0.0001 ~ "<0.0001",
                                 p.adj >= 0.0001 & p.adj < 0.001 ~ as.character(round(p.adj, 4)),
                                 p.adj >= 0.001 & p.adj < 0.01  ~ as.character(round(p.adj, 3)),
                                 p.adj >= 0.01              ~ as.character(round(p.adj, 2)),
                                 .default = NA)) %>% 
                               rename(Candidate = Metabolite,
                                      Step2_p_adj = p.adj),
                             by = "Candidate")

## Beurre

### Moyenne (95% CI)
tableau_beurre <- left_join(x = tableau_beurre,
                            y = data_beurre_delta %>%
                              select(all_of(metabolite_beurre)) %>%
                              get_summary_stats(type = "mean_ci") %>%
                              mutate(mean = round(mean, 2),
                                     ci = round(ci, 2),
                                     mean_ci = paste0(mean, " (", (mean - ci), ", ", (mean + ci), ")")) %>%
                              select(variable, mean_ci) %>%
                              rename(Candidate = variable,
                                     Step2_mean_ci = mean_ci),
                            by = "Candidate")

### p-value
tableau_beurre <- left_join(x = tableau_beurre,
                            y = test_sign_beurre %>%
                              select(Metabolite, p.adj) %>%
                              mutate(p.adj = case_when(
                                p.adj < 0.0001 ~ "<0.0001",
                                p.adj >= 0.0001 & p.adj < 0.001 ~ as.character(round(p.adj, 4)),
                                p.adj >= 0.001 & p.adj < 0.01  ~ as.character(round(p.adj, 3)),
                                p.adj >= 0.01              ~ as.character(round(p.adj, 2)),
                                .default = NA)) %>% 
                              rename(Candidate = Metabolite,
                                     Step2_p_adj = p.adj),
                            by = "Candidate")


# Extraction des métabolites significatifs (passent step 2)
metabolite_lait <- subset(test_sign_lait, subset = p.adj < 0.05)$Metabolite
metabolite_fromage <- subset(test_sign_fromage, subset = p.adj < 0.05)$Metabolite
metabolite_beurre <- subset(test_sign_beurre, subset = p.adj < 0.05)$Metabolite

rm(test_sign_lait, test_sign_fromage, test_sign_beurre)

# Step 3 - Test statistique de la specificité ----

# Fonction
test_spec <- function(biomarkers, data_long, data_delta, intensity_col = "Intensity",
                      treatment_col = "Traitement", metabolite_col = "Metabolite", p_threshold = 0.05) {
  
  # Test de Wilcoxon des biomarqueurs de l'aliment cible dans les données de l'aliment interférent
  test_res <- data_long %>%
    filter(!!sym(metabolite_col) %in% biomarkers) %>%
    group_by(!!sym(metabolite_col)) %>%
    wilcox_test(as.formula(paste(intensity_col, "~", treatment_col)), paired = TRUE, alternative = "less") %>%
    adjust_pvalue(method = "none") %>%
    add_significance()
  
  # Calcul des médianes dans data_delta (assume : colonnes = métabolites, lignes = sujets ou paires)
  med_df <- data.frame(Mediane = apply(data_delta[, biomarkers, drop = FALSE], 2, median)) %>%
    rownames_to_column(var = metabolite_col)
  
  # Fusion des résultats statistiques et des médianes
  test_res <- merge(test_res, med_df, by = metabolite_col)
  
  # Retourner uniquement les métabolites non significatifs ou de médiane ≤ 0
  biomarkers_filtered <- test_res %>%
    filter(!(p < p_threshold)) %>%
    pull(!!sym(metabolite_col))
  
  biomarkers_unspecific <- test_res %>%
    filter(p.adj < p_threshold) %>%
    pull(!!sym(metabolite_col))
  
  
  return(list(
    resultats_complets = test_res,
    biomarqueurs_specifiques = biomarkers_filtered,
    biomarkers_non_specifiques = biomarkers_unspecific
  ))
}

# Lait
test_spec_lait_fromage <- test_spec(metabolite_lait, data_fromage_long, data_fromage_delta)
test_spec_lait_beurre <- test_spec(metabolite_lait, data_beurre_long, data_beurre_delta)

# Fromage
test_spec_fromage_lait <- test_spec(metabolite_fromage, data_lait_long, data_lait_delta)
test_spec_fromage_beurre <- test_spec(metabolite_fromage, data_beurre_long, data_beurre_delta)

# Beurre
test_spec_beurre_lait <- test_spec(metabolite_beurre, data_lait_long, data_lait_delta)
test_spec_beurre_fromage <- test_spec(metabolite_beurre, data_fromage_long, data_fromage_delta)


# Extraction des résultats pour tableau supplémentaire

## Lait

### Moyenne (95% CI)

# Specificité fromage #
tableau_lait <- left_join(x = tableau_lait,
                          y = data_fromage_delta %>%
                            select(all_of(metabolite_lait)) %>%
                            get_summary_stats(type = "mean_ci") %>%
                            mutate(mean = round(mean, 2),
                                   ci = round(ci, 2),
                                   mean_ci = paste0(mean, " (", (mean - ci), ", ", (mean + ci), ")")) %>%
                            select(variable, mean_ci) %>%
                            rename(Candidate = variable,
                                   Step3_cheese_mean_ci = mean_ci),
                          by = "Candidate")

# Specificité beurre #
tableau_lait <- left_join(x = tableau_lait,
                          y = data_beurre_delta %>%
                            select(all_of(metabolite_lait)) %>%
                            get_summary_stats(type = "mean_ci") %>%
                            mutate(mean = round(mean, 2),
                                   ci = round(ci, 2),
                                   mean_ci = paste0(mean, " (", (mean - ci), ", ", (mean + ci), ")")) %>%
                            select(variable, mean_ci) %>%
                            rename(Candidate = variable,
                                   Step3_butter_mean_ci = mean_ci),
                          by = "Candidate")

### p-value

# Specificité fromage #
tableau_lait <- left_join(x = tableau_lait,
                          y = test_spec_lait_fromage$resultats_complets %>%
                            select(Metabolite, p.adj) %>%
                            mutate(p.adj = case_when(
                              p.adj < 0.0001 ~ "<0.0001",
                              p.adj >= 0.0001 & p.adj < 0.001 ~ as.character(round(p.adj, 4)),
                              p.adj >= 0.001 & p.adj < 0.01  ~ as.character(round(p.adj, 3)),
                              p.adj >= 0.01              ~ as.character(round(p.adj, 2)),
                              .default = NA)) %>% 
                            rename(Candidate = Metabolite,
                                   Step3_cheese_p_adj = p.adj),
                          by = "Candidate")

# Specificité beurre #
tableau_lait <- left_join(x = tableau_lait,
                          y = test_spec_lait_beurre$resultats_complets %>%
                            select(Metabolite, p.adj) %>%
                            mutate(p.adj = case_when(
                              p.adj < 0.0001 ~ "<0.0001",
                              p.adj >= 0.0001 & p.adj < 0.001 ~ as.character(round(p.adj, 4)),
                              p.adj >= 0.001 & p.adj < 0.01  ~ as.character(round(p.adj, 3)),
                              p.adj >= 0.01              ~ as.character(round(p.adj, 2)),
                              .default = NA)) %>% 
                            rename(Candidate = Metabolite,
                                   Step3_butter_p_adj = p.adj),
                          by = "Candidate")

## Fromage

### Moyenne (95% CI)

# Specificité lait #
tableau_fromage <- left_join(x = tableau_fromage,
                             y = data_lait_delta %>%
                               select(all_of(metabolite_fromage)) %>%
                               get_summary_stats(type = "mean_ci") %>%
                               mutate(mean = round(mean, 2),
                                      ci = round(ci, 2),
                                      mean_ci = paste0(mean, " (", (mean - ci), ", ", (mean + ci), ")")) %>%
                               select(variable, mean_ci) %>%
                               rename(Candidate = variable,
                                      Step3_milk_mean_ci = mean_ci),
                             by = "Candidate")

# Specificité beurre #
tableau_fromage <- left_join(x = tableau_fromage,
                             y = data_beurre_delta %>%
                               select(all_of(metabolite_fromage)) %>%
                               get_summary_stats(type = "mean_ci") %>%
                               mutate(mean = round(mean, 2),
                                      ci = round(ci, 2),
                                      mean_ci = paste0(mean, " (", (mean - ci), ", ", (mean + ci), ")")) %>%
                               select(variable, mean_ci) %>%
                               rename(Candidate = variable,
                                      Step3_butter_mean_ci = mean_ci),
                             by = "Candidate")

### p-value

# Specificité lait #
tableau_fromage <- left_join(x = tableau_fromage,
                             y = test_spec_fromage_lait$resultats_complets %>%
                               select(Metabolite, p.adj) %>%
                               mutate(p.adj = case_when(
                                 p.adj < 0.0001 ~ "<0.0001",
                                 p.adj >= 0.0001 & p.adj < 0.001 ~ as.character(round(p.adj, 4)),
                                 p.adj >= 0.001 & p.adj < 0.01  ~ as.character(round(p.adj, 3)),
                                 p.adj >= 0.01              ~ as.character(round(p.adj, 2)),
                                 .default = NA)) %>% 
                               rename(Candidate = Metabolite,
                                      Step3_milk_p_adj = p.adj),
                             by = "Candidate")

# Specificité beurre #
tableau_fromage <- left_join(x = tableau_fromage,
                             y = test_spec_fromage_beurre$resultats_complets %>%
                               select(Metabolite, p.adj) %>%
                               mutate(p.adj = case_when(
                                 p.adj < 0.0001 ~ "<0.0001",
                                 p.adj >= 0.0001 & p.adj < 0.001 ~ as.character(round(p.adj, 4)),
                                 p.adj >= 0.001 & p.adj < 0.01  ~ as.character(round(p.adj, 3)),
                                 p.adj >= 0.01              ~ as.character(round(p.adj, 2)),
                                 .default = NA)) %>% 
                               rename(Candidate = Metabolite,
                                      Step3_butter_p_adj = p.adj),
                             by = "Candidate")



## Beurre

### Moyenne (95% CI)

# Specificité fromage #
tableau_beurre <- left_join(x = tableau_beurre,
                            y = data_fromage_delta %>%
                              select(all_of(metabolite_beurre)) %>%
                              get_summary_stats(type = "mean_ci") %>%
                              mutate(mean = round(mean, 2),
                                     ci = round(ci, 2),
                                     mean_ci = paste0(mean, " (", (mean - ci), ", ", (mean + ci), ")")) %>%
                              select(variable, mean_ci) %>%
                              rename(Candidate = variable,
                                     Step3_cheese_mean_ci = mean_ci),
                            by = "Candidate")

# Specificité lait #
tableau_beurre <- left_join(x = tableau_beurre,
                            y = data_lait_delta %>%
                              select(all_of(metabolite_beurre)) %>%
                              get_summary_stats(type = "mean_ci") %>%
                              mutate(mean = round(mean, 2),
                                     ci = round(ci, 2),
                                     mean_ci = paste0(mean, " (", (mean - ci), ", ", (mean + ci), ")")) %>%
                              select(variable, mean_ci) %>%
                              rename(Candidate = variable,
                                     Step3_milk_mean_ci = mean_ci),
                            by = "Candidate")
### p-value

# Specificité fromage #
tableau_beurre <- left_join(x = tableau_beurre,
                            y = test_spec_beurre_fromage$resultats_complets %>%
                              select(Metabolite, p.adj) %>%
                              mutate(p.adj = case_when(
                                p.adj < 0.0001 ~ "<0.0001",
                                p.adj >= 0.0001 & p.adj < 0.001 ~ as.character(round(p.adj, 4)),
                                p.adj >= 0.001 & p.adj < 0.01  ~ as.character(round(p.adj, 3)),
                                p.adj >= 0.01              ~ as.character(round(p.adj, 2)),
                                .default = NA)) %>% 
                              rename(Candidate = Metabolite,
                                     Step3_cheese_p_adj = p.adj),
                            by = "Candidate")

# Specificité lait #
tableau_beurre <- left_join(x = tableau_beurre,
                            y = test_spec_beurre_lait$resultats_complets %>%
                              select(Metabolite, p.adj) %>%
                              mutate(p.adj = case_when(
                                p.adj < 0.0001 ~ "<0.0001",
                                p.adj >= 0.0001 & p.adj < 0.001 ~ as.character(round(p.adj, 4)),
                                p.adj >= 0.001 & p.adj < 0.01  ~ as.character(round(p.adj, 3)),
                                p.adj >= 0.01              ~ as.character(round(p.adj, 2)),
                                .default = NA)) %>% 
                              rename(Candidate = Metabolite,
                                     Step3_milk_p_adj = p.adj),
                            by = "Candidate")



# Extraction des métabolites

## Lait
metabolite_lait <- intersect(test_spec_lait_fromage$biomarqueurs_specifiques,
                             test_spec_lait_beurre$biomarqueurs_specifiques)

## Fromage
metabolite_fromage <- intersect(test_spec_fromage_lait$biomarqueurs_specifiques,
                                test_spec_fromage_beurre$biomarqueurs_specifiques)

## Beurre
metabolite_beurre <- intersect(test_spec_beurre_lait$biomarqueurs_specifiques,
                               test_spec_beurre_fromage$biomarqueurs_specifiques)

## Total dairy
metabolite_total_dairy <- unique(c(
  intersect(test_spec_lait_fromage$biomarkers_non_specifiques, test_spec_lait_beurre$biomarkers_non_specifiques), # Lait
  intersect(test_spec_fromage_lait$biomarkers_non_specifiques, test_spec_fromage_beurre$biomarkers_non_specifiques), # Fromage
  intersect(test_spec_beurre_lait$biomarkers_non_specifiques, test_spec_beurre_fromage$biomarkers_non_specifiques))) # Beurre

## Dairy fat
metabolite_dairy_fat <- setdiff(unique(c(setdiff(test_spec_fromage_beurre$biomarkers_non_specifiques, test_spec_fromage_lait$biomarkers_non_specifiques), # Fromage
                                 setdiff(test_spec_beurre_fromage$biomarkers_non_specifiques, test_spec_beurre_lait$biomarkers_non_specifiques))), # Beurre
                                metabolite_total_dairy)

# Tableau supplémentaire dairy fat
tableau_dairy_fat <- full_join(x = full_join(x = data_fromage_delta %>%
                                               select(all_of(metabolite_dairy_fat)) %>%
                                               get_summary_stats(type = "mean_ci") %>%
                                               mutate(mean = round(mean, 2),
                                                      ci = round(ci, 2),
                                                      mean_ci = paste0(mean, " (", (mean - ci), ", ", (mean + ci), ")")) %>%
                                               select(variable, mean_ci) %>%
                                               rename(Candidate = variable,
                                                      Step3_cheese_mean_ci = mean_ci),
                                             y = subset(data_fromage_long, subset = Metabolite %in% metabolite_dairy_fat) %>%
                                               group_by(Metabolite) %>%
                                               wilcox_test(Intensity ~ Traitement, paired = TRUE, alternative = "less") %>%
                                               adjust_pvalue(method = "none") %>%
                                               select(Metabolite, p.adj) %>%
                                               mutate(p.adj = case_when(
                                                 p.adj < 0.0001 ~ "<0.0001",
                                                 p.adj >= 0.0001 & p.adj < 0.001 ~ as.character(round(p.adj, 4)),
                                                 p.adj >= 0.001 & p.adj < 0.01  ~ as.character(round(p.adj, 3)),
                                                 p.adj >= 0.01              ~ as.character(round(p.adj, 2)),
                                                 .default = NA)) %>%
                                               rename("Candidate" = "Metabolite",
                                                      "Step3_cheese_p_adj" = "p.adj"),
                                             by = "Candidate"),
                               y = full_join(x = data_beurre_delta %>%
                                               select(all_of(metabolite_dairy_fat)) %>%
                                               get_summary_stats(type = "mean_ci") %>%
                                               mutate(mean = round(mean, 2),
                                                      ci = round(ci, 2),
                                                      mean_ci = paste0(mean, " (", (mean - ci), ", ", (mean + ci), ")")) %>%
                                               select(variable, mean_ci) %>%
                                               rename(Candidate = variable,
                                                      Step3_butter_mean_ci = mean_ci),
                                             y = subset(data_beurre_long, subset = Metabolite %in% metabolite_dairy_fat) %>%
                                               group_by(Metabolite) %>%
                                               wilcox_test(Intensity ~ Traitement, paired = TRUE, alternative = "less") %>%
                                               adjust_pvalue(method = "none") %>%
                                               select(Metabolite, p.adj) %>%
                                               mutate(p.adj = case_when(
                                                 p.adj < 0.0001 ~ "<0.0001",
                                                 p.adj >= 0.0001 & p.adj < 0.001 ~ as.character(round(p.adj, 4)),
                                                 p.adj >= 0.001 & p.adj < 0.01  ~ as.character(round(p.adj, 3)),
                                                 p.adj >= 0.01              ~ as.character(round(p.adj, 2)),
                                                 .default = NA)) %>%
                                               rename("Candidate" = "Metabolite",
                                                      "Step3_butter_p_adj" = "p.adj"),
                                             by = "Candidate"),
                               by = "Candidate")



# Step 4 - Test statistique pour la robustesse en interventionnelle - Projet PLI ----

# Lait
test_robust_lait <- subset(data_pli_long, subset = Metabolite %in% metabolite_lait) %>%
  group_by(Metabolite) %>%
  wilcox_test(Intensity ~ Traitement, paired = TRUE, alternative = "less") %>%
  adjust_pvalue(method = "none") %>%
  add_significance()

# Fromage
test_robust_fromage <- subset(data_pli_long, subset = Metabolite %in% metabolite_fromage) %>%
  group_by(Metabolite) %>%
  wilcox_test(Intensity ~ Traitement, paired = TRUE, alternative = "less") %>%
  adjust_pvalue(method = "none") %>%
  add_significance()

# Dairy fat
test_robust_dairy_fat <- subset(data_pli_long, subset = Metabolite %in% metabolite_dairy_fat) %>%
  group_by(Metabolite) %>%
  wilcox_test(Intensity ~ Traitement, paired = TRUE, alternative = "less") %>%
  adjust_pvalue(method = "none") %>%
  add_significance()


# Extraction des résultats pour tableau supplémentaire

## Lait

### Moyenne (95% CI)
tableau_lait <- left_join(x = tableau_lait,
                          y = data_pli_delta %>%
                            select(all_of(metabolite_lait)) %>%
                            get_summary_stats(type = "mean_ci") %>%
                            mutate(mean = round(mean, 2),
                                   ci = round(ci, 2),
                                   mean_ci = paste0(mean, " (", (mean - ci), ", ", (mean + ci), ")")) %>%
                            select(variable, mean_ci) %>%
                            rename(Candidate = variable,
                                   Step4_mean_ci = mean_ci),
                          by = "Candidate")

### p-value
tableau_lait <- left_join(x = tableau_lait,
                          y = test_robust_lait %>%
                            select(Metabolite, p.adj) %>%
                            mutate(p.adj = case_when(
                              p.adj < 0.0001 ~ "<0.0001",
                              p.adj >= 0.0001 & p.adj < 0.001 ~ as.character(round(p.adj, 4)),
                              p.adj >= 0.001 & p.adj < 0.01  ~ as.character(round(p.adj, 3)),
                              p.adj >= 0.01              ~ as.character(round(p.adj, 2)),
                              .default = NA)) %>% 
                            rename(Candidate = Metabolite,
                                   Step4_p_adj = p.adj),
                          by = "Candidate")

## Fromage

### Moyenne (95% CI)
tableau_fromage <- left_join(x = tableau_fromage,
                             y = data_pli_delta %>%
                               select(all_of(metabolite_fromage)) %>%
                               get_summary_stats(type = "mean_ci") %>%
                               mutate(mean = round(mean, 2),
                                      ci = round(ci, 2),
                                      mean_ci = paste0(mean, " (", (mean - ci), ", ", (mean + ci), ")")) %>%
                               select(variable, mean_ci) %>%
                               rename(Candidate = variable,
                                      Step4_mean_ci = mean_ci),
                             by = "Candidate")

### p-value
tableau_fromage <- left_join(x = tableau_fromage,
                             y = test_robust_fromage %>%
                               select(Metabolite, p.adj) %>%
                               mutate(p.adj = case_when(
                                 p.adj < 0.0001 ~ "<0.0001",
                                 p.adj >= 0.0001 & p.adj < 0.001 ~ as.character(round(p.adj, 4)),
                                 p.adj >= 0.001 & p.adj < 0.01  ~ as.character(round(p.adj, 3)),
                                 p.adj >= 0.01              ~ as.character(round(p.adj, 2)),
                                 .default = NA)) %>% 
                               rename(Candidate = Metabolite,
                                      Step4_p_adj = p.adj),
                             by = "Candidate")

## Dairy fat

### Moyenne (95% CI)
tableau_dairy_fat <- left_join(x = tableau_dairy_fat,
                               y = data_pli_delta %>%
                                 select(all_of(metabolite_dairy_fat)) %>%
                                 get_summary_stats(type = "mean_ci") %>%
                                 mutate(mean = round(mean, 2),
                                        ci = round(ci, 2),
                                        mean_ci = paste0(mean, " (", (mean - ci), ", ", (mean + ci), ")")) %>%
                                 select(variable, mean_ci) %>%
                                 rename(Candidate = variable,
                                        Step4_mean_ci = mean_ci),
                               by = "Candidate")

### p-value
tableau_dairy_fat <- left_join(x = tableau_dairy_fat,
                               y = test_robust_dairy_fat %>%
                                 select(Metabolite, p.adj) %>%
                                 mutate(p.adj = case_when(
                                   p.adj < 0.0001 ~ "<0.0001",
                                   p.adj >= 0.0001 & p.adj < 0.001 ~ as.character(round(p.adj, 4)),
                                   p.adj >= 0.001 & p.adj < 0.01  ~ as.character(round(p.adj, 3)),
                                   p.adj >= 0.01              ~ as.character(round(p.adj, 2)),
                                   .default = NA)) %>% 
                                 rename(Candidate = Metabolite,
                                        Step4_p_adj = p.adj),
                               by = "Candidate")



# Extraction des métabolites significatifs
metabolite_lait <- subset(test_robust_lait, subset = p.adj < 0.05)$Metabolite
metabolite_fromage <- subset(test_robust_fromage, subset = p.adj < 0.05)$Metabolite
metabolite_dairy_fat <- subset(test_robust_dairy_fat, subset = p.adj < 0.05)$Metabolite

# Sauvegarde des features sélectionnés ----
save(metabolite_lait, metabolite_fromage, metabolite_beurre, metabolite_dairy_fat,
     file = "Metabo/Features Post Step 4.rda")

# Sauvegarde des tableaux supplémentaires ----
save(tableau_lait, tableau_fromage, tableau_beurre, tableau_dairy_fat,
     file = "Metabo/Tableaux S7-S8-S9-S10 - Steps 2-3-4.rda")
