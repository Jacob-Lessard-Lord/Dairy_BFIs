# Date : 16 octobre 2025
# Auteur : Jacob Lessard-Lord

# Projet PLC

# Analyse dose-response RCT

# Untargeted 

# Set working directory ----
setwd("C:/Users/jackl/OneDrive/Bureau/École/Postdoctorat/Chercheur/Jean-Philippe/PLC/Processing_PLC/Untargeted/Final_Processing/Data")

# Chargement des packages requis ----
library(tidyverse)
library(ggpubr)
library(rstatix)
library(emmeans)
library(openxlsx)

# Importation des données métabolomiques ----
load("Metabo/Untargeted - ALL.rda")

# Importation  des biomarqueurs ----
load("Metabo/Features Post Step 4.rda")

# Mise en forme des données ----

# Création d'un vecteur avec tous les biomarqueurs
all_biomarkers <- c(metabolite_lait, metabolite_fromage, metabolite_beurre, metabolite_dairy_fat)

# Selection des données post-dietes RCT
data_ech <- subset(data_ech, subset = Traitement %in% c("Lait", "Fromage", "Beurre", "Controle", "Produits laitiers"),
                   select = c("Projet", "Sujet", "Projet_Sujet", "Traitement", all_biomarkers))

# Scaling
data_ech[,  all_biomarkers] <- scale(log(data_ech[,  all_biomarkers]))

# Creation des datasets ----

# Lait
data_lait <- rbind(data_ech %>% filter(Projet == "HDL" & Traitement == "Controle") %>% mutate(Lait_groupe = "Q1"),
                   data_ech %>% filter(Projet == "PLI" & Traitement == "Produits laitiers") %>% mutate(Lait_groupe = "Q2"),
                   data_ech %>% filter(Traitement == "Lait") %>% mutate(Lait_groupe = "Q3")) %>%
  mutate(Lait_groupe = ordered(Lait_groupe))

# Fromage
data_fromage <- rbind(data_ech %>% filter(Projet == "LAIT" & Traitement == "Controle") %>% mutate(Fromage_groupe = "Q1"),
                      data_ech %>% filter(Projet == "PLI" & Traitement == "Produits laitiers") %>% mutate(Fromage_groupe = "Q2"),
                      data_ech %>% filter(Traitement == "Fromage") %>% mutate(Fromage_groupe = "Q3")) %>%
  mutate(Fromage_groupe = ordered(Fromage_groupe))

# Dairy fat
data_dairy_fat <- rbind(data_ech %>% filter(Projet == "LAIT" & Traitement == "Controle") %>% mutate(Dairy_fat_groupe = "Q1"),
                        data_ech %>% filter(Projet == "PLI" & Traitement == "Produits laitiers") %>% mutate(Dairy_fat_groupe = "Q2"),
                        data_ech %>% filter(Traitement %in% c("Fromage", "Beurre")) %>% mutate(Dairy_fat_groupe = "Q3")) %>%
  mutate(Dairy_fat_groupe = ordered(Dairy_fat_groupe)) %>% 
  group_by(Projet_Sujet, Dairy_fat_groupe) %>%
  summarise(across(any_of(all_biomarkers), ~ mean(.x, na.rm = TRUE)), .groups = "drop")


# Fonction pour régression du facteur ordonné ----
reg_fun <- function(data_df, metabolite, cat_aliment) {
  
  # Entrainement du modele de regression
  modele <- lm(as.formula(paste("Intensity ~", cat_aliment)), data = subset(data_df, subset = Metabolite == metabolite))
  
  # Extraction des resultats de la regression
  results_reg <- summary(modele)
  
  # Calcul de la moyenne pour chacun des quartiles
  quartiles_moy <- as.data.frame(emmeans(modele, reformulate(cat_aliment))) %>%
    mutate(mean_95ci = paste0(round(emmean, 2), " (", round(lower.CL, 2), ", ", round(upper.CL, 2), ")")) # Mise en forme de la moyenne (95% CI)
  
  # Enregistrement des resultats
  results_df <- data.frame(Metabolite = metabolite,
                           Coef_beta = results_reg$coefficients[2, "Estimate"],
                           p_value = results_reg$coefficients[2, "Pr(>|t|)"],
                           r2 = results_reg$r.squared,
                           Q1 = quartiles_moy[1, ]$emmean,
                           Q2 = quartiles_moy[2, ]$emmean,
                           Q3 = quartiles_moy[3, ]$emmean,
                           Q1_95ci = quartiles_moy[1, ]$mean_95ci,
                           Q2_95ci = quartiles_moy[2, ]$mean_95ci,
                           Q3_95ci = quartiles_moy[3, ]$mean_95ci)
  
  return(results_df)
  
}

# Fonction pour tester le dose-response ----
cohort_criteria_fun <- function(data_df,
                                metabolites,
                                group_col,
                                tag,
                                rename_map = NULL,
                                alpha = 0.05,
                                do_plot = TRUE) {
  
  # data         : data.frame (format large)
  # metabolites  : vecteur des colonnes métabolites à tester
  # group_col    : nom de la variable de groupe (ex. "Lait_groupe")
  # tag          : suffixe pour nommer les colonnes (ex. "ffq", "r24w_q3", ...)
  # rename_map   : nommage lisible des métabolites (named vector: old="New name")
  # alpha        : seuil p-value ajusté BH
  # do_plot      : TRUE/FALSE pour retourner un ggplot
  
  # Stockage du nom de la variable de groupe sous forme d'expression à évaluer
  grp <- rlang::sym(group_col)
  
  # Long
  long <- data_df %>%
    select(all_of(c(group_col, unname(metabolites)))) %>%
    pivot_longer(cols = all_of(unname(metabolites)), names_to = "Metabolite", values_to = "Intensity")
  
  
  # Régression par quartiles sur les métabolites retenus
  reg_df <- if (length(metabolites)) {
    purrr::map_dfr(metabolites, ~ reg_fun(long, .x, group_col))
  } else {
    tibble()
  }
  
  if (nrow(reg_df)) {
    reg_df <- reg_df %>% mutate(p_value_adj = p.adjust(p_value, method = "none"))
    # Critère strict : Q1 < Q2 < Q3
    met_final <- reg_df %>%
      filter(p_value_adj < alpha, Q1 < Q2, Q2 < Q3) %>%
      pull(Metabolite)
  } else {
    reg_df <- tibble(Metabolite = character(), p_value = numeric(), p_value_adj = numeric())
    met_final <- character()
  }
  
  reg_out <- reg_df %>%
    select(Metabolite, p_value_adj) %>%
    rename(!!rlang::sym(paste0(tag, "_groupe")) := p_value_adj)
  
  # Graphique (optionnel) sur les métabolites finaux
  plt <- NULL
  if (do_plot && length(met_final)) {
    plot_data <- long %>% filter(Metabolite %in% met_final)
    if (!is.null(rename_map)) {
      plot_data <- plot_data %>% mutate(Metabolite = recode(Metabolite, !!!rename_map))
    }
    plt <- ggerrorplot(data = plot_data,
                       x = "Metabolite", y = "Intensity", color = group_col,
                       desc_stat = "mean_ci", ylab = "Scaled area", xlab = "Biomarker",
                       palette = c("#01903d", "#507bff", "#8a335d")) +
      #scale_x_discrete(labels = setNames(names(metabolites), unname(metabolites))) +
      coord_flip() +
      theme(legend.title = element_blank())
  }
  
  list(
    results = reg_out,
    biomarqueurs = met_final,      # après régression + monotonie stricte
    plot = plt
  )
}

# Analyse dose-response ----

# Biomarqueurs lait
results_lait <- cohort_criteria_fun(data_df = data_lait,
                                    metabolites = metabolite_lait,
                                    group_col = "Lait_groupe",
                                    tag = "Trial")

# Biomarqueurs fromage
results_fromage <- cohort_criteria_fun(data_df = data_fromage,
                                       metabolites = metabolite_fromage,
                                       group_col = "Fromage_groupe",
                                       tag = "Trial")

# Biomarqueurs dairy fat
results_dairy_fat <- cohort_criteria_fun(data_df = data_dairy_fat,
                                         metabolites = metabolite_dairy_fat,
                                         group_col = "Dairy_fat_groupe",
                                         tag = "Trial")


# Mise en forme des résultats pour tableaux supplémentaires ----

# Lait
tableau_lait <- full_join(x = data_lait %>%
                            select(all_of(c("Lait_groupe", metabolite_lait))) %>%
                            pivot_longer(cols = all_of(metabolite_lait), names_to = "Candidate", values_to = "Intensity") %>%
                            group_by(Candidate, Lait_groupe) %>%
                            get_summary_stats(Intensity, type = "mean_ci") %>%
                            mutate(mean = round(mean, 2),
                                   ci = round(ci, 2),
                                   mean_ci = paste0(mean, " (", (mean - ci), ", ", (mean + ci), ")")) %>%
                            select(Lait_groupe, Candidate, mean_ci) %>%
                            pivot_wider(names_from = Lait_groupe, values_from = mean_ci) %>%
                            rename("Step5_D1" = "Q1",
                                   "Step5_D2" = "Q2",
                                   "Step5_D3" = "Q3"),
                          y = results_lait$results %>% rename("Candidate" = "Metabolite", "Step5_p_value" = "Trial_groupe"),
                          by = "Candidate") %>%
  mutate(Step5_p_value = case_when(Step5_p_value < 0.0001 ~ "<0.0001",
                                  Step5_p_value >= 0.0001 & Step5_p_value < 0.001 ~ as.character(round(Step5_p_value, 4)),
                                  Step5_p_value >= 0.001 & Step5_p_value < 0.01  ~ as.character(round(Step5_p_value, 3)),
                                  Step5_p_value >= 0.01 ~ as.character(round(Step5_p_value, 2)),
                                  .default = NA))

    
# Fromage
tableau_fromage <- full_join(x = data_fromage %>%
                               select(all_of(c("Fromage_groupe", metabolite_fromage))) %>%
                               pivot_longer(cols = all_of(metabolite_fromage), names_to = "Candidate", values_to = "Intensity") %>%
                               group_by(Candidate, Fromage_groupe) %>%
                               get_summary_stats(Intensity, type = "mean_ci") %>%
                               mutate(mean = round(mean, 2),
                                      ci = round(ci, 2),
                                      mean_ci = paste0(mean, " (", (mean - ci), ", ", (mean + ci), ")")) %>%
                               select(Fromage_groupe, Candidate, mean_ci) %>%
                               pivot_wider(names_from = Fromage_groupe, values_from = mean_ci) %>%
                               rename("Step5_D1" = "Q1",
                                      "Step5_D2" = "Q2",
                                      "Step5_D3" = "Q3"),
                             y = results_fromage$results %>% rename("Candidate" = "Metabolite", "Step5_p_value" = "Trial_groupe"),
                             by = "Candidate") %>%
  mutate(Step5_p_value = case_when(Step5_p_value < 0.0001 ~ "<0.0001",
                                   Step5_p_value >= 0.0001 & Step5_p_value < 0.001 ~ as.character(round(Step5_p_value, 4)),
                                   Step5_p_value >= 0.001 & Step5_p_value < 0.01  ~ as.character(round(Step5_p_value, 3)),
                                   Step5_p_value >= 0.01 ~ as.character(round(Step5_p_value, 2)),
                                   .default = NA))

# Dairy fat
tableau_dairy_fat <- full_join(x = data_dairy_fat %>%
                                 select(all_of(c("Dairy_fat_groupe", metabolite_dairy_fat))) %>%
                                 pivot_longer(cols = all_of(metabolite_dairy_fat), names_to = "Candidate", values_to = "Intensity") %>%
                                 group_by(Candidate, Dairy_fat_groupe) %>%
                                 get_summary_stats(Intensity, type = "mean_ci") %>%
                                 mutate(mean = round(mean, 2),
                                        ci = round(ci, 2),
                                        mean_ci = paste0(mean, " (", (mean - ci), ", ", (mean + ci), ")")) %>%
                                 select(Dairy_fat_groupe, Candidate, mean_ci) %>%
                                 pivot_wider(names_from = Dairy_fat_groupe, values_from = mean_ci) %>%
                                 rename("Step5_D1" = "Q1",
                                        "Step5_D2" = "Q2",
                                        "Step5_D3" = "Q3"),
                               y = results_dairy_fat$results %>% rename("Candidate" = "Metabolite", "Step5_p_value" = "Trial_groupe"),
                               by = "Candidate") %>%
  mutate(Step5_p_value = case_when(Step5_p_value < 0.0001 ~ "<0.0001",
                                   Step5_p_value >= 0.0001 & Step5_p_value < 0.001 ~ as.character(round(Step5_p_value, 4)),
                                   Step5_p_value >= 0.001 & Step5_p_value < 0.01  ~ as.character(round(Step5_p_value, 3)),
                                   Step5_p_value >= 0.01 ~ as.character(round(Step5_p_value, 2)),
                                   .default = NA))



# Extraction des biomarqueurs ----
metabolite_lait <- results_lait$biomarqueurs
metabolite_fromage <- results_fromage$biomarqueurs
metabolite_dairy_fat <- results_dairy_fat$biomarqueurs

# Sauvegarde des features sélectionnés ----
save(metabolite_lait, metabolite_fromage, metabolite_beurre, metabolite_dairy_fat,
     file = "Metabo/Features Post Step 5.rda")

# Sauvegarde des tableaux supplémentaires ----
save(tableau_lait, tableau_fromage, tableau_dairy_fat,
     file = "Metabo/Tableaux S7-S8-S10 - Step 5.rda")
