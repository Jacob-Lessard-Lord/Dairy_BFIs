# Date : 19 novembre 2025
# Auteur : Jacob Lessard-Lord

# Projet PLC

# ÉValuation de la performance des biomarqueurs en contexte observationnel

# Niveau de dVB plasmatique en fonction de l'apport auto-rapporté de lait et de yogourt

# Untargeted 

# Set working directory ----
setwd("C:/Users/jackl/OneDrive/Bureau/École/Postdoctorat/Chercheur/Jean-Philippe/PLC/Processing_PLC/Untargeted/Final_Processing/Data")

# Chargement des packages requis ----
library(tidyverse)
library(ggpubr)
library(rstatix)
library(emmeans)

# Importation des données nutritionnelles ----
load("FFQ/Data FFQ.rda")
load("R24W/Data R24W.rda")

# Importation des données métabolomiques ----
load("Metabo/Untargeted - ALL.rda")

# Importation  des biomarqueurs ----
load("Metabo/Biomarqueurs finaux.rda")

# Mise en forme des données ----

# Création d'un vecteur avec tous les biomarqueurs
all_biomarkers <- c(biomarqueur_lait, biomarqueur_beurre, biomarqueur_dairy_fat)

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



# Analyse statistique ----
results_yogourt_dvb <- cohort_criteria_fun(data_df = data_r24w_mean %>% filter(Lait_groupe == "Q1 - Very low consumption (<15g)"),
                                           metabolites = "RPLC_Pos_0.75_160.1302mz", # dVB
                                           group_col = "Yogourt_groupe",
                                           tag = "r24w_mean")


# Graphique ----
graph_r24w_mean_yogourt <- ggboxplot(data = data_r24w_mean %>% filter(Lait_groupe == "Q1 - Very low consumption (<15g)") %>%
                                       mutate(Yogourt_groupe = factor(case_when(Yogourt == 0 ~ "0 g/day",
                                                                                Yogourt > 0 & Yogourt <= 125 ~ "0-125 g/day",
                                                                                Yogourt > 125 ~ ">125 g/day"),
                                                                      levels = c("0 g/day", "0-125 g/day", ">125 g/day"))),
                                     x = "Yogourt_groupe", y = "RPLC_Pos_0.75_160.1302mz", add = "jitter", fill = "Yogourt_groupe",
                                     palette = c("#01903d", "#507bff", "#8a335d"),
                                     xlab = "Self-reported yogurt consumption", ylab = "δ-Valerobetaine\nnormalized abundance") +
  theme(legend.position = "none")

# Exportation de la figure ----
setwd("C:/Users/jackl/OneDrive/Bureau/École/Postdoctorat/Chercheur/Jean-Philippe/PLC/Article/Biomarqueurs/Figures")

tiff(filename = "Figure S3 - dVB Yogourt.tiff", width = 10, height = 12, units = "cm", res = 600, compression = "lzw")

graph_r24w_mean_yogourt

dev.off()


