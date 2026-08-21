# Date : 13 juillet 2026
# Auteur : Jacob Lessard-Lord

# Projet PLC

# ÉValuation de la performance des biomarqueurs en contexte observationnel

# Reproductibilité à court et long term

# Untargeted 

# Set working directory ----
setwd("C:/Users/jackl/OneDrive/Bureau/École/Postdoctorat/Chercheur/Jean-Philippe/PLC/Processing_PLC/Untargeted/Final_Processing/Data")

# Chargement des packages requis ----
library(tidyverse)
library(ggpubr)
library(rstatix)
library(openxlsx)
library(irr)

# Importation des données ----
load("FH/Data FFQ FH.rda")

# Mise en forme des données ----

# Calcul temps entre visite FH-LPA-CAC
data_ffq_fh_cac_lpa %>% 
  select(Projet, nomat, date_visite) %>%
  mutate(nomat = case_when(Projet == "FH-CAC" ~ str_remove(nomat, "^1"),
                           Projet == "FH-LPA" ~ nomat),
         Projet = str_replace(Projet, "-", "_")) %>%
  pivot_wider(names_from = "Projet", values_from = "date_visite") %>%
  mutate(Time_between_visit = as.numeric(FH_CAC - FH_LPA) / 365.25) %>%
  get_summary_stats(Time_between_visit, type = "common")

# Mise en forme longue
data_ffq_fh_diet <- data_ffq_fh_diet %>%
  select(nomat, Visites,
         delta_valerobetaine, Homostachydrine, PE_20_4, PC_20_4, SM_35_1, PC_35_1, SM_43_2, SM_43_1, Cer_43_1,
         Lait, Beurre, Dairy_fat) %>%
  pivot_longer(cols = c(delta_valerobetaine, Homostachydrine, PE_20_4, PC_20_4, SM_35_1, PC_35_1, SM_43_2, SM_43_1, Cer_43_1,
                        Lait, Beurre, Dairy_fat),
               names_to = "Exposition", values_to = "Valeur") %>%
  pivot_wider(names_from = "Visites", values_from = "Valeur")

data_ffq_fh_cac_lpa <- data_ffq_fh_cac_lpa %>%
  select(Projet, nomat,
         delta_valerobetaine, Homostachydrine, PE_20_4, PC_20_4, SM_35_1, PC_35_1, SM_43_2, SM_43_1, Cer_43_1,
         Lait, Beurre, Dairy_fat) %>%
  mutate(nomat = case_when(Projet == "FH-CAC" ~ str_remove(nomat, "^1"),
                           Projet == "FH-LPA" ~ nomat),
         Projet = str_replace(Projet, "-", "_")) %>%
  pivot_longer(cols = c(delta_valerobetaine, Homostachydrine, PE_20_4, PC_20_4, SM_35_1, PC_35_1, SM_43_2, SM_43_1, Cer_43_1,
                        Lait, Beurre, Dairy_fat),
               names_to = "Exposition", values_to = "Valeur") %>%
  pivot_wider(names_from = "Projet", values_from = "Valeur")

# Calcul des ICC ----

# Reproductibilité à court terme
icc_fh_diet <- data_ffq_fh_diet %>%
  group_by(Exposition) %>%
  group_split() %>%
  map_dfr(function(sub_df) {
    
    # Extraction du nom de l'exposition
    exp_name <- unique(sub_df$Exposition)
    
    # Mise en forme en matrice
    mat <- sub_df %>% select(V0, V1) %>% as.data.frame()
    
    # Conserver seulement les paires complètes
    mat <- mat[complete.cases(mat), ]
    
    # Calcul de l'ICC
    res <- icc(mat, model = "twoway", type = "agreement", unit = "single")
    
    # Mise en forme des résultats
    tibble(
      Exposition = exp_name,
      n = res$subjects,
      icc = res$value,
      lower_ci = res$lbound,
      upper_ci = res$ubound,
      p_value = res$p.value
    )
  })

# Reproductibilité à long terme
icc_fh_cac_lpa <- data_ffq_fh_cac_lpa %>%
  group_by(Exposition) %>%
  group_split() %>%
  map_dfr(function(sub_df) {
    
    # Extraction du nom de l'exposition
    exp_name <- unique(sub_df$Exposition)
    
    # Mise en forme en matrice
    mat <- sub_df %>% select(FH_LPA, FH_CAC) %>% as.data.frame()
    
    # Conserver seulement les paires complètes
    mat <- mat[complete.cases(mat), ]
    
    # Calcul de l'ICC
    res <- icc(mat, model = "twoway", type = "agreement", unit = "single")
    
    # Mise en forme des résultats
    tibble(
      Exposition = exp_name,
      n = res$subjects,
      icc = res$value,
      lower_ci = res$lbound,
      upper_ci = res$ubound,
      p_value = res$p.value
    )
  })


# Calcul des RDR (Regression dilution ratio) ----

# Reproductibilité à court terme 
rdr_fh_diet <- data_ffq_fh_diet %>%
  group_by(Exposition) %>%
  group_split() %>%
  map_dfr(function(sub_df) {
    
    # Extraction du nom de l'exposition
    exp_name <- unique(sub_df$Exposition)
    
    # Mise en forme en matrice
    mat <- sub_df %>% select(V0, V1) %>% as.data.frame()
    mat <- mat[complete.cases(mat), ]
    
    # Nombre de paires
    n_pairs <- nrow(mat)
    
    # Regression dilution ratio = pente de V1 ~ V0
    fit <- lm(V1 ~ V0, data = mat)
    slope <- coef(fit)["V0"]
    ci <- confint(fit)["V0", ]
    pval <- summary(fit)$coefficients["V0", "Pr(>|t|)"]
    r2 <- summary(fit)$r.squared
    
    # Mise en forme des résultats
    tibble(
      Exposition = exp_name,
      n = n_pairs,
      rdr = slope,
      lower_ci = ci[1],
      upper_ci = ci[2],
      p_value = pval,
      r_squared = r2
    )
  })

# Reproductibilité à long terme 
rdr_fh_cac_lpa <- data_ffq_fh_cac_lpa %>%
  group_by(Exposition) %>%
  group_split() %>%
  map_dfr(function(sub_df) {
    
    # Extraction du nom de l'exposition
    exp_name <- unique(sub_df$Exposition)
    
    # Mise en forme en matrice
    mat <- sub_df %>% select(FH_LPA, FH_CAC) %>% as.data.frame()
    mat <- mat[complete.cases(mat), ]
    
    # Nombre de paires
    n_pairs <- nrow(mat)
    
    # Regression dilution ratio = pente de FH_CAC (2e mesure) ~ FH_LPA (1ere mesure)
    fit <- lm(FH_CAC ~ FH_LPA, data = mat)
    slope <- coef(fit)["FH_LPA"]
    ci <- confint(fit)["FH_LPA", ]
    pval <- summary(fit)$coefficients["FH_LPA", "Pr(>|t|)"]
    r2 <- summary(fit)$r.squared
    
    # Mise en forme des résultats
    tibble(
      Exposition = exp_name,
      n = n_pairs,
      rdr = slope,
      lower_ci = ci[1],
      upper_ci = ci[2],
      p_value = pval,
      r_squared = r2
    )
  })

# Mise en forme finale ----
data_icc_rdr <- full_join(x = full_join(x = icc_fh_diet %>%
                                          mutate(icc = sprintf("%.2f", icc),
                                                 lower_ci = sprintf("%.2f", lower_ci),
                                                 upper_ci = sprintf("%.2f", upper_ci),
                                                 icc_short_term = paste0(icc, " (", lower_ci, ", ", upper_ci, ")")) %>%
                                          select(Exposition, icc_short_term),
                                        y = rdr_fh_diet %>%
                                          mutate(rdr = sprintf("%.2f", rdr),
                                                 lower_ci = sprintf("%.2f", lower_ci),
                                                 upper_ci = sprintf("%.2f", upper_ci),
                                                 rdr_short_term = paste0(rdr, " (", lower_ci, ", ", upper_ci, ")")) %>%
                                          select(Exposition, rdr_short_term),
                                        by = "Exposition"),
                          y = full_join(x = icc_fh_cac_lpa %>%
                                          mutate(icc = sprintf("%.2f", icc),
                                                 lower_ci = sprintf("%.2f", lower_ci),
                                                 upper_ci = sprintf("%.2f", upper_ci),
                                                 icc_long_term = paste0(icc, " (", lower_ci, ", ", upper_ci, ")")) %>%
                                          select(Exposition, icc_long_term),
                                        y = rdr_fh_cac_lpa %>%
                                          mutate(rdr = sprintf("%.2f", rdr),
                                                 lower_ci = sprintf("%.2f", lower_ci),
                                                 upper_ci = sprintf("%.2f", upper_ci),
                                                 rdr_long_term = paste0(rdr, " (", lower_ci, ", ", upper_ci, ")")) %>%
                                          select(Exposition, rdr_long_term),
                                        by = "Exposition"),
                          by = "Exposition") %>%
  mutate(Food = case_when(Exposition %in% c("Lait", "delta_valerobetaine", "Homostachydrine") ~ "Milk",
                          Exposition %in% c("Beurre", "PE_20_4", "PC_20_4") ~ "Butter",
                          Exposition %in% c("Dairy_fat", "SM_35_1", "PC_35_1", "SM_43_2", "SM_43_1", "Cer_43_1") ~ "Concentrated full-fat dairy"),
         Exposition = case_match(Exposition,
                                 "Beurre" ~ "FFQ",
                                 "Cer_43_1" ~ "Cer d18:1/25:0",
                                 "Dairy_fat" ~ "FFQ",
                                 "Homostachydrine" ~ "Homostachydrine",
                                 "Lait" ~ "FFQ",
                                 "PC_20_4" ~ "PC (P-17:0/20:4, O-17:1/20:4)",
                                 "PC_35_1" ~ "PC 17:0/18:1",
                                 "PE_20_4" ~ "PE (P-17:0/20:4, O-17:1/20:4)",
                                 "SM_35_1" ~ "SM d18:1/17:0",
                                 "SM_43_1" ~ "SM 43:1",
                                 "SM_43_2" ~ "SM 43:2",
                                 "delta_valerobetaine" ~ "δ-Valerobetaine",
                                 .default = Exposition))

# Exportation des résultats ----
write.xlsx(data_icc_rdr, file = "Metabo/Tables SX - Short-term & long-term reproducibility.xlsx")










