# Date : 24 juillet 2025
# Auteur : Jacob Lessard-Lord

# Projet PLC

# ÉValuation de la performance des biomarqueurs en contexte observationnel

# Wilcox test / Régression ordinale

# Untargeted 

# Set working directory ----
setwd("C:/Users/jackl/OneDrive/Bureau/École/Postdoctorat/Chercheur/Jean-Philippe/PLC/Processing_PLC/Untargeted/Final_Processing/Data")

# Chargement des packages requis ----
library(tidyverse)
library(ggpubr)
library(rstatix)
library(emmeans)
library(readxl)
library(openxlsx)

# Importation des données nutritionnelles ----
load("FFQ/Data FFQ.rda")
load("R24W/Data R24W.rda")

# Importation des données métabolomiques ----
load("Metabo/Untargeted - ALL.rda")

# Importation  des biomarqueurs ----
load("Metabo/Biomarqueurs post-dereplication.rda")

# Mise en forme des données ----

# Création d'un vecteur avec tous les biomarqueurs
all_biomarkers <- c(metabolite_lait, metabolite_fromage, metabolite_beurre, metabolite_dairy_fat)

# Selection des données au baseline des RCT + eMECA
data_ech <- subset(data_ech, Traitement == "Baseline" | Projet == "eMECA", select = c("Projet", "Sujet", "Projet_Sujet", all_biomarkers))

# Scaling
data_ech[,  all_biomarkers] <- scale(log(data_ech[,  all_biomarkers]))

# Merge avec les données de FFQ
data_ffq <- merge(x = data_ffq, y = data_ech, by = c("Projet", "Sujet"))

# Merge avec les données de R24W
data_r24w_q3 <- merge(x = data_r24w_q3, y = data_ech, by = c("Projet", "Sujet"))

data_r24w_mean <- merge(x = data_r24w_mean, y = data_ech, by = c("Projet", "Sujet"))

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

# Initialisation des tableaux de résultats pour l'article ----

# Fonction pour tester trend linéaire Q1 vs Q2 vs Q3 ----
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
    rename(!!rlang::sym(paste0(tag, "_quartile")) := p_value_adj)
  
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
      #coord_flip() +
      theme(legend.title = element_blank())
  }
  
  list(
    results = reg_out,
    biomarqueurs = met_final,      # après régression + monotonie stricte
    plot = plt
  )
}

# FFQ ----

# Lait
lait_ffq <- cohort_criteria_fun(data_df = data_ffq,
                                metabolites = metabolite_lait,
                                group_col = "Lait_groupe",
                                tag = "ffq")

# Fromage
# n = 0

# Beurre
beurre_ffq <- cohort_criteria_fun(data_df = data_ffq,
                                   metabolites = metabolite_beurre,
                                   group_col = "Beurre_groupe",
                                   tag = "ffq")

# Dairy fat
dairy_fat_ffq <- cohort_criteria_fun(data_df = data_ffq,
                                  metabolites = metabolite_dairy_fat,
                                  group_col = "Dairy_fat_groupe",
                                  tag = "ffq")

# R24W - Moyenne des 3 questionnaires - Tous les sujets ----

# Lait
lait_r24w_mean <- cohort_criteria_fun(data_df = data_r24w_mean,
                                      metabolites = metabolite_lait,
                                      group_col = "Lait_groupe",
                                      tag = "r24w_mean")

# Fromage
# n = 0

# Beurre
beurre_r24w_mean <- cohort_criteria_fun(data_df = data_r24w_mean,
                                        metabolites = metabolite_beurre,
                                        group_col = "Beurre_groupe",
                                        tag = "r24w_mean")

# Dairy fat
dairy_fat_r24w_mean <- cohort_criteria_fun(data_df = data_r24w_mean,
                                           metabolites = metabolite_dairy_fat,
                                           group_col = "Dairy_fat_groupe",
                                           tag = "r24w_mean")

# R24W - Q3 ----

# Lait
lait_r24w_q3 <- cohort_criteria_fun(data_df = data_r24w_q3,
                                    metabolites = metabolite_lait,
                                    group_col = "Lait_groupe",
                                    tag = "r24w_q3")

# Fromage
# n = 0

# Beurre
beurre_r24w_q3 <- cohort_criteria_fun(data_df = data_r24w_q3,
                                      metabolites = metabolite_beurre,
                                      group_col = "Beurre_groupe",
                                      tag = "r24w_q3")

# Dairy fat
dairy_fat_r24w_q3 <- cohort_criteria_fun(data_df = data_r24w_q3,
                                         metabolites = metabolite_dairy_fat,
                                         group_col = "Dairy_fat_groupe",
                                         tag = "r24w_q3")


# Mise en forme des résultats ----

#### Lait ####

# Merge des résultats par aliment
tableau_lait <- full_join(x = full_join(x = lait_ffq$results,
                                        y = lait_r24w_mean$results,
                                        by = "Metabolite"),
                          y = lait_r24w_q3$results,
                          by = "Metabolite")

# Mise en forme des p-values
tableau_lait <- tableau_lait %>%
  mutate(across(
    .cols = -Metabolite,
    .fns = ~ case_when(
      . < 0.0001 ~ "<0.0001",
      . >= 0.0001 & . < 0.001 ~ as.character(round(., 4)),
      . >= 0.001 & . < 0.01  ~ as.character(round(., 3)),
      . >= 0.01              ~ as.character(round(., 2)),
      TRUE                   ~ NA_character_  # pour conserver les NA
    )
  ))

# Calcul des moyennes (95% CI) par quartile

## FFQ

### Initialisation du data.frame
results_lait_ffq <- data.frame(Metabolite = NA, Coef_beta = NA, p_value = NA, r2 = NA,
                               Q1 = NA, Q2 = NA, Q3 = NA, 
                               Q1_95ci = NA, Q2_95ci = NA, Q3_95ci = NA)

### Calcul
for(i in 1:length(metabolite_lait)) {
  results_lait_ffq[i, ] <- reg_fun(data_ffq %>%
                                     select(all_of(c("Lait_groupe", unname(metabolite_lait)))) %>%
                                     pivot_longer(cols = all_of(unname(metabolite_lait)), names_to = "Metabolite", values_to = "Intensity"),
                                   metabolite_lait[i], "Lait_groupe")
}

### Mise en forme
results_lait_ffq <- results_lait_ffq %>%
  select(Metabolite, Q1_95ci, Q2_95ci, Q3_95ci) %>%
  rename(Q1_ffq = Q1_95ci, Q2_ffq = Q2_95ci, Q3_ffq = Q3_95ci)

## R24W - Moyenne des 3 questionnaires

### Initialisation du data.frame
results_lait_r24w_mean <- data.frame(Metabolite = NA, Coef_beta = NA, p_value = NA, r2 = NA,
                                     Q1 = NA, Q2 = NA, Q3 = NA, 
                                     Q1_95ci = NA, Q2_95ci = NA, Q3_95ci = NA)

### Calcul
for(i in 1:length(metabolite_lait)) {
  results_lait_r24w_mean[i, ] <- reg_fun(data_r24w_mean %>%
                                           select(all_of(c("Lait_groupe", unname(metabolite_lait)))) %>%
                                           pivot_longer(cols = all_of(unname(metabolite_lait)), names_to = "Metabolite", values_to = "Intensity"),
                                         metabolite_lait[i], "Lait_groupe")
}

### Mise en forme
results_lait_r24w_mean <- results_lait_r24w_mean %>%
  select(Metabolite, Q1_95ci, Q2_95ci, Q3_95ci) %>%
  rename(Q1_r24w_mean = Q1_95ci, Q2_r24w_mean = Q2_95ci, Q3_r24w_mean = Q3_95ci)

## R24 - Questionnaire #3

### Initialisation du data.frame
results_lait_r24w_q3 <- data.frame(Metabolite = NA, Coef_beta = NA, p_value = NA, r2 = NA,
                               Q1 = NA, Q2 = NA, Q3 = NA, 
                               Q1_95ci = NA, Q2_95ci = NA, Q3_95ci = NA)

### Calcul
for(i in 1:length(metabolite_lait)) {
  results_lait_r24w_q3[i, ] <- reg_fun(data_r24w_q3 %>%
                                         select(all_of(c("Lait_groupe", unname(metabolite_lait)))) %>%
                                         pivot_longer(cols = all_of(unname(metabolite_lait)), names_to = "Metabolite", values_to = "Intensity"),
                                       metabolite_lait[i], "Lait_groupe")
}

### Mise en forme
results_lait_r24w_q3 <- results_lait_r24w_q3 %>%
  select(Metabolite, Q1_95ci, Q2_95ci, Q3_95ci) %>%
  rename(Q1_r24w_q3 = Q1_95ci, Q2_r24w_q3 = Q2_95ci, Q3_r24w_q3 = Q3_95ci)

# Merge de tous les résultats
tableau_lait <- full_join(x = full_join(x = full_join(x = full_join(x = data.frame(ID = names(metabolite_lait),
                                                                                                               Metabolite = metabolite_lait),
                                                                                                y = tableau_lait,
                                                                                                by = "Metabolite"),
                                                                                  y = results_lait_ffq,
                                                                                  by = "Metabolite"),
                                                                    y = results_lait_r24w_mean,
                                                                    by = "Metabolite"),
                                                      y = results_lait_r24w_q3,
                                                      by = "Metabolite") %>%
  select(ID, Metabolite, Q1_ffq, Q2_ffq, Q3_ffq, ffq_quartile,
         Q1_r24w_mean, Q2_r24w_mean, Q3_r24w_mean, r24w_mean_quartile,
         Q1_r24w_q3, Q2_r24w_q3, Q3_r24w_q3, r24w_q3_quartile)

rm(results_lait_ffq, results_lait_r24w_mean, results_lait_r24w_q3)


#### Fromage ####
# n = 0

#### Beurre ####

# Merge des résultats par aliment
tableau_beurre <- full_join(x = full_join(x = beurre_ffq$results,
                                        y = beurre_r24w_mean$results,
                                        by = "Metabolite"),
                          y = beurre_r24w_q3$results,
                          by = "Metabolite")

# Mise en forme des p-values
tableau_beurre <- tableau_beurre %>%
  mutate(across(
    .cols = -Metabolite,
    .fns = ~ case_when(
      . < 0.0001 ~ "<0.0001",
      . >= 0.0001 & . < 0.001 ~ as.character(round(., 4)),
      . >= 0.001 & . < 0.01  ~ as.character(round(., 3)),
      . >= 0.01              ~ as.character(round(., 2)),
      TRUE                   ~ NA_character_  # pour conserver les NA
    )
  ))

# Calcul des moyennes (95% CI) par quartile

## FFQ

### Initialisation du data.frame
results_beurre_ffq <- data.frame(Metabolite = NA, Coef_beta = NA, p_value = NA, r2 = NA,
                               Q1 = NA, Q2 = NA, Q3 = NA, 
                               Q1_95ci = NA, Q2_95ci = NA, Q3_95ci = NA)

### Calcul
for(i in 1:length(metabolite_beurre)) {
  results_beurre_ffq[i, ] <- reg_fun(data_ffq %>%
                                     select(all_of(c("Beurre_groupe", unname(metabolite_beurre)))) %>%
                                     pivot_longer(cols = all_of(unname(metabolite_beurre)), names_to = "Metabolite", values_to = "Intensity"),
                                   metabolite_beurre[i], "Beurre_groupe")
}

### Mise en forme
results_beurre_ffq <- results_beurre_ffq %>%
  select(Metabolite, Q1_95ci, Q2_95ci, Q3_95ci) %>%
  rename(Q1_ffq = Q1_95ci, Q2_ffq = Q2_95ci, Q3_ffq = Q3_95ci)

## R24W - Moyenne des 3 questionnaires

### Initialisation du data.frame
results_beurre_r24w_mean <- data.frame(Metabolite = NA, Coef_beta = NA, p_value = NA, r2 = NA,
                                     Q1 = NA, Q2 = NA, Q3 = NA, 
                                     Q1_95ci = NA, Q2_95ci = NA, Q3_95ci = NA)

### Calcul
for(i in 1:length(metabolite_beurre)) {
  results_beurre_r24w_mean[i, ] <- reg_fun(data_r24w_mean %>%
                                           select(all_of(c("Beurre_groupe", unname(metabolite_beurre)))) %>%
                                           pivot_longer(cols = all_of(unname(metabolite_beurre)), names_to = "Metabolite", values_to = "Intensity"),
                                         metabolite_beurre[i], "Beurre_groupe")
}

### Mise en forme
results_beurre_r24w_mean <- results_beurre_r24w_mean %>%
  select(Metabolite, Q1_95ci, Q2_95ci, Q3_95ci) %>%
  rename(Q1_r24w_mean = Q1_95ci, Q2_r24w_mean = Q2_95ci, Q3_r24w_mean = Q3_95ci)

## R24 - Questionnaire #3

### Initialisation du data.frame
results_beurre_r24w_q3 <- data.frame(Metabolite = NA, Coef_beta = NA, p_value = NA, r2 = NA,
                                   Q1 = NA, Q2 = NA, Q3 = NA, 
                                   Q1_95ci = NA, Q2_95ci = NA, Q3_95ci = NA)

### Calcul
for(i in 1:length(metabolite_beurre)) {
  results_beurre_r24w_q3[i, ] <- reg_fun(data_r24w_q3 %>%
                                         select(all_of(c("Beurre_groupe", unname(metabolite_beurre)))) %>%
                                         pivot_longer(cols = all_of(unname(metabolite_beurre)), names_to = "Metabolite", values_to = "Intensity"),
                                       metabolite_beurre[i], "Beurre_groupe")
}

### Mise en forme
results_beurre_r24w_q3 <- results_beurre_r24w_q3 %>%
  select(Metabolite, Q1_95ci, Q2_95ci, Q3_95ci) %>%
  rename(Q1_r24w_q3 = Q1_95ci, Q2_r24w_q3 = Q2_95ci, Q3_r24w_q3 = Q3_95ci)

# Merge de tous les résultats
tableau_beurre <- full_join(x = full_join(x = full_join(x = full_join(x = data.frame(ID = names(metabolite_beurre),
                                                                                   Metabolite = metabolite_beurre),
                                                                    y = tableau_beurre,
                                                                    by = "Metabolite"),
                                                      y = results_beurre_ffq,
                                                      by = "Metabolite"),
                                        y = results_beurre_r24w_mean,
                                        by = "Metabolite"),
                          y = results_beurre_r24w_q3,
                          by = "Metabolite") %>%
  select(ID, Metabolite, Q1_ffq, Q2_ffq, Q3_ffq, ffq_quartile,
         Q1_r24w_mean, Q2_r24w_mean, Q3_r24w_mean, r24w_mean_quartile,
         Q1_r24w_q3, Q2_r24w_q3, Q3_r24w_q3, r24w_q3_quartile)

rm(results_beurre_ffq, results_beurre_r24w_mean, results_beurre_r24w_q3)

#### Dairy fat ####

# Merge des résultats par aliment
tableau_dairy_fat <- full_join(x = full_join(x = dairy_fat_ffq$results,
                                        y = dairy_fat_r24w_mean$results,
                                        by = "Metabolite"),
                          y = dairy_fat_r24w_q3$results,
                          by = "Metabolite")

# Mise en forme des p-values
tableau_dairy_fat <- tableau_dairy_fat %>%
  mutate(across(
    .cols = -Metabolite,
    .fns = ~ case_when(
      . < 0.0001 ~ "<0.0001",
      . >= 0.0001 & . < 0.001 ~ as.character(round(., 4)),
      . >= 0.001 & . < 0.01  ~ as.character(round(., 3)),
      . >= 0.01              ~ as.character(round(., 2)),
      TRUE                   ~ NA_character_  # pour conserver les NA
    )
  ))

# Calcul des moyennes (95% CI) par quartile

## FFQ

### Initialisation du data.frame
results_dairy_fat_ffq <- data.frame(Metabolite = NA, Coef_beta = NA, p_value = NA, r2 = NA,
                               Q1 = NA, Q2 = NA, Q3 = NA, 
                               Q1_95ci = NA, Q2_95ci = NA, Q3_95ci = NA)

### Calcul
for(i in 1:length(metabolite_dairy_fat)) {
  results_dairy_fat_ffq[i, ] <- reg_fun(data_ffq %>%
                                     select(all_of(c("Dairy_fat_groupe", unname(metabolite_dairy_fat)))) %>%
                                     pivot_longer(cols = all_of(unname(metabolite_dairy_fat)), names_to = "Metabolite", values_to = "Intensity"),
                                   metabolite_dairy_fat[i], "Dairy_fat_groupe")
}

### Mise en forme
results_dairy_fat_ffq <- results_dairy_fat_ffq %>%
  select(Metabolite, Q1_95ci, Q2_95ci, Q3_95ci) %>%
  rename(Q1_ffq = Q1_95ci, Q2_ffq = Q2_95ci, Q3_ffq = Q3_95ci)

## R24W - Moyenne des 3 questionnaires

### Initialisation du data.frame
results_dairy_fat_r24w_mean <- data.frame(Metabolite = NA, Coef_beta = NA, p_value = NA, r2 = NA,
                                     Q1 = NA, Q2 = NA, Q3 = NA, 
                                     Q1_95ci = NA, Q2_95ci = NA, Q3_95ci = NA)

### Calcul
for(i in 1:length(metabolite_dairy_fat)) {
  results_dairy_fat_r24w_mean[i, ] <- reg_fun(data_r24w_mean %>%
                                           select(all_of(c("Dairy_fat_groupe", unname(metabolite_dairy_fat)))) %>%
                                           pivot_longer(cols = all_of(unname(metabolite_dairy_fat)), names_to = "Metabolite", values_to = "Intensity"),
                                         metabolite_dairy_fat[i], "Dairy_fat_groupe")
}

### Mise en forme
results_dairy_fat_r24w_mean <- results_dairy_fat_r24w_mean %>%
  select(Metabolite, Q1_95ci, Q2_95ci, Q3_95ci) %>%
  rename(Q1_r24w_mean = Q1_95ci, Q2_r24w_mean = Q2_95ci, Q3_r24w_mean = Q3_95ci)

## R24 - Questionnaire #3

### Initialisation du data.frame
results_dairy_fat_r24w_q3 <- data.frame(Metabolite = NA, Coef_beta = NA, p_value = NA, r2 = NA,
                                   Q1 = NA, Q2 = NA, Q3 = NA, 
                                   Q1_95ci = NA, Q2_95ci = NA, Q3_95ci = NA)

### Calcul
for(i in 1:length(metabolite_dairy_fat)) {
  results_dairy_fat_r24w_q3[i, ] <- reg_fun(data_r24w_q3 %>%
                                         select(all_of(c("Dairy_fat_groupe", unname(metabolite_dairy_fat)))) %>%
                                         pivot_longer(cols = all_of(unname(metabolite_dairy_fat)), names_to = "Metabolite", values_to = "Intensity"),
                                       metabolite_dairy_fat[i], "Dairy_fat_groupe")
}

### Mise en forme
results_dairy_fat_r24w_q3 <- results_dairy_fat_r24w_q3 %>%
  select(Metabolite, Q1_95ci, Q2_95ci, Q3_95ci) %>%
  rename(Q1_r24w_q3 = Q1_95ci, Q2_r24w_q3 = Q2_95ci, Q3_r24w_q3 = Q3_95ci)

# Merge de tous les résultats
tableau_dairy_fat <- full_join(x = full_join(x = full_join(x = full_join(x = data.frame(ID = names(metabolite_dairy_fat),
                                                                                   Metabolite = metabolite_dairy_fat),
                                                                    y = tableau_dairy_fat,
                                                                    by = "Metabolite"),
                                                      y = results_dairy_fat_ffq,
                                                      by = "Metabolite"),
                                        y = results_dairy_fat_r24w_mean,
                                        by = "Metabolite"),
                          y = results_dairy_fat_r24w_q3,
                          by = "Metabolite") %>%
  select(ID, Metabolite, Q1_ffq, Q2_ffq, Q3_ffq, ffq_quartile,
         Q1_r24w_mean, Q2_r24w_mean, Q3_r24w_mean, r24w_mean_quartile,
         Q1_r24w_q3, Q2_r24w_q3, Q3_r24w_q3, r24w_q3_quartile)

rm(results_dairy_fat_ffq, results_dairy_fat_r24w_mean, results_dairy_fat_r24w_q3)

# Exportation des résultats ----
tableau_lait <- tableau_lait %>% rename("Candidate" = "Metabolite",
                                        "Step6_D1_ffq" = "Q1_ffq", "Step6_D2_ffq" = "Q2_ffq", "Step6_D3_ffq" = "Q3_ffq", "Step_6_p_value_ffq" = "ffq_quartile",
                                        "Step6_D1_r24w_mean" = "Q1_r24w_mean", "Step6_D2_r24w_mean" = "Q2_r24w_mean", "Step6_D3_r24w_mean" = "Q3_r24w_mean",
                                        "Step_6_p_value_r24w_mean" = "r24w_mean_quartile",
                                        "Step6_D1_r24w_q3" = "Q1_r24w_q3", "Step6_D2_r24w_q3" = "Q2_r24w_q3", "Step6_D3_r24w_q3" = "Q3_r24w_q3",
                                        "Step_6_p_value_r24w_q3" = "r24w_q3_quartile")

tableau_beurre <- tableau_beurre %>% rename("Candidate" = "Metabolite",
                                            "Step6_D1_ffq" = "Q1_ffq", "Step6_D2_ffq" = "Q2_ffq", "Step6_D3_ffq" = "Q3_ffq", "Step_6_p_value_ffq" = "ffq_quartile",
                                            "Step6_D1_r24w_mean" = "Q1_r24w_mean", "Step6_D2_r24w_mean" = "Q2_r24w_mean", "Step6_D3_r24w_mean" = "Q3_r24w_mean",
                                            "Step_6_p_value_r24w_mean" = "r24w_mean_quartile",
                                            "Step6_D1_r24w_q3" = "Q1_r24w_q3", "Step6_D2_r24w_q3" = "Q2_r24w_q3", "Step6_D3_r24w_q3" = "Q3_r24w_q3",
                                            "Step_6_p_value_r24w_q3" = "r24w_q3_quartile")

tableau_dairy_fat <- tableau_dairy_fat %>% rename("Candidate" = "Metabolite",
                                                  "Step6_D1_ffq" = "Q1_ffq", "Step6_D2_ffq" = "Q2_ffq", "Step6_D3_ffq" = "Q3_ffq", "Step_6_p_value_ffq" = "ffq_quartile",
                                                  "Step6_D1_r24w_mean" = "Q1_r24w_mean", "Step6_D2_r24w_mean" = "Q2_r24w_mean", "Step6_D3_r24w_mean" = "Q3_r24w_mean",
                                                  "Step_6_p_value_r24w_mean" = "r24w_mean_quartile",
                                                  "Step6_D1_r24w_q3" = "Q1_r24w_q3", "Step6_D2_r24w_q3" = "Q2_r24w_q3", "Step6_D3_r24w_q3" = "Q3_r24w_q3",
                                                  "Step_6_p_value_r24w_q3" = "r24w_q3_quartile")

save(tableau_lait, tableau_beurre,tableau_dairy_fat,
     file = "Metabo/Tableaux S7-S9-S10 - Step 6.rda")

# Enregistrement des metabolites ----

# Lait
biomarqueur_lait <- unique(c(lait_ffq$biomarqueurs, lait_r24w_mean$biomarqueurs))

# Fromage
# n = 0

# Beurre
biomarqueur_beurre <- unique(c(beurre_ffq$biomarqueurs, beurre_r24w_mean$biomarqueurs))

# Dairy fat
biomarqueur_dairy_fat <- unique(c(dairy_fat_ffq$biomarqueurs, dairy_fat_r24w_mean$biomarqueurs))

# Enregistrement
save(biomarqueur_lait, biomarqueur_beurre, biomarqueur_dairy_fat,
     file = "Metabo/Biomarqueurs finaux.rda")

# Enregistrer les résultats ----
save(lait_ffq, lait_r24w_mean, lait_r24w_q3,
     beurre_ffq, beurre_r24w_mean, beurre_r24w_q3,
     dairy_fat_ffq, dairy_fat_r24w_mean, dairy_fat_r24w_q3,
     file = "Metabo/Results Step 6.rda")
