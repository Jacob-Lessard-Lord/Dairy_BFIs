# Date : 15 août 2025
# Auteur : Jacob Lessard-Lord

# Projet PLC

# Données sur les facteurs de rsique

# Pre-processing

# Chargement des packages requis ----
library(tidyverse)
library(rstatix)
library(readxl)
library(openxlsx)

# Set working directory ----
setwd("C:/Users/jackl/OneDrive/Bureau/École/Postdoctorat/Chercheur/Jean-Philippe/PLC/Processing_PLC/Untargeted/Final_Processing/Data")

# Importation des données ----

# Données de FANI
data_fani <- read_excel("FANI/2025-07-22 - PLC - Rapport FANI.xlsx", guess_max = 10000, skip = 3, na = ".")

# Données de metabo
load("Metabo/Untargeted - ALL.rda")

# Données de FFQ
load("FFQ/Data FFQ.rda")
load("R24W/Data R24W.rda")

# Mise en forme des données ----

# Sélection des données
data_clinique <- rbind(
  
  # LAIT 
  data_fani %>%
    filter(Projet == "LAIT" & Visites == "s0") %>%
    select(c(Projet, nomat, age, sex, bmi)),
  
  # GABA2
  data_fani %>%
    filter(Projet == "GABA2-" & Visites == "Scr 2") %>%
    select(c(Projet, nomat, age, sexe_ncp, bmi)) %>% 
    rename(sex = sexe_ncp) %>%
    mutate(Projet = str_replace(Projet, "GABA2-", "GABA2")),
  
  # HDL
  data_fani %>%
    filter(Projet == "HDL" & Visites == "Scr 1") %>%
    select(c(Projet, nomat, age, sex, bmi)),
  
  # PLI
  data_fani %>%
    filter(Projet == "PLI" & Visites == "S1") %>%
    select(c(Projet, nomat, age, sexe_ncp, bmi)) %>% 
    rename(sex = sexe_ncp),
  
  # eMECA
  full_join(full_join(data_fani %>%
                        filter(Projet == "eMECA-" & Visites == "S1") %>%
                        select(c(Projet, nomat, sexe_ncp)) %>% 
                        rename(sex = sexe_ncp), 
                      data_fani %>%
                        filter(Projet == "eMECA-" & Visites == "V1") %>%
                        select(c(Projet, nomat, bmi)),
                      by = c("Projet", "nomat")), 
            data_fani %>%
              filter(Projet == "eMECA-" & Visites == "Questweb1") %>%
              select(c(Projet, nomat, age)),
            by = c("Projet", "nomat")) %>%
    mutate(Projet = str_replace(Projet, "eMECA-", "eMECA")) %>%
    select(Projet, nomat, age, sex, bmi)) %>% 
  mutate(sex = factor(case_match(sex, 1 ~ "Female", 2 ~ "Male")),
         Sujet = factor(as.numeric(nomat)),
         Projet_Sujet = paste0(Projet, "_", Sujet)) %>%
  filter(Projet_Sujet %in% unique(data_ech$Projet_Sujet))

# Mise en forme données de FFQ et R24W ----

# FFQ
data_ffq <- data_ffq %>%
  mutate(Sujet = factor(Sujet),
         Projet_Sujet = paste0(Projet, "_", Sujet))

# R24W
data_r24w_mean <- data_r24w_mean %>%
  mutate(Sujet = factor(Sujet),
         Projet_Sujet = paste0(Projet, "_", Sujet))

# Ajout de données de FFQ ----
data_clinique <- rbind(
  
  # LAIT
  merge(x = data_clinique %>% filter(Projet == "LAIT"),
        y = data_ffq %>% select(Projet_Sujet, Lait, Fromage, Beurre, Dairy_fat, Total_dairy),
        by = "Projet_Sujet", all.x = TRUE),
  
  # GABA2
  merge(x = data_clinique %>% filter(Projet == "GABA2"),
        y = data_ffq %>% select(Projet_Sujet, Lait, Fromage, Beurre, Dairy_fat, Total_dairy),
        by = "Projet_Sujet", all.x = TRUE),
  
  # HDL
  merge(x = data_clinique %>% filter(Projet == "HDL"),
        y = data_ffq %>% select(Projet_Sujet, Lait, Fromage, Beurre, Dairy_fat, Total_dairy),
        by = "Projet_Sujet", all.x = TRUE),
  
  # PLI 
  merge(x = data_clinique %>% filter(Projet == "PLI"),
        y = data_ffq %>% select(Projet_Sujet, Lait, Fromage, Beurre, Dairy_fat, Total_dairy),
        by = "Projet_Sujet", all.x = TRUE),
  
  # eMECA
  merge(x = merge(x = data_clinique %>% filter(Projet == "eMECA"),
                  y = data_r24w_mean %>% select(Projet_Sujet),
                  by = "Projet_Sujet"),
        y = data_ffq %>% select(Projet_Sujet, Lait, Fromage, Beurre, Dairy_fat, Total_dairy),
        by = "Projet_Sujet", all.x = TRUE))

# n par trial ----

# LAIT
nrow(data_clinique %>% filter(Projet == "LAIT"))

# GABA2 
nrow(data_clinique %>% filter(Projet == "GABA2"))

# HDL
nrow(data_clinique %>% filter(Projet == "HDL"))

# PLI
nrow(data_clinique %>% filter(Projet == "PLI"))

# eMECA
nrow(data_clinique %>% filter(Projet == "eMECA"))

# FFQ baseline
nrow(data_clinique %>% filter(Projet != "eMECA") %>% filter(!is.na(Lait)))

# Calcul des statistiques descriptives ----
results <- rbind(data_clinique %>%
  group_by(Projet) %>%
  get_summary_stats(type = "mean_sd") %>%
  mutate(Results_all = paste(round(mean, 1), " ± ", round(sd, 1))) %>%
  select(-c(n, mean, sd)),
  data_clinique %>% freq_table(Projet, sex) %>% mutate(Results_all = paste0(n, " (", prop, "%)")) %>% select(-c(n, prop)) %>% rename (variable = sex)) %>%
  pivot_wider(names_from = "Projet", values_from = "Results_all") %>%
  select(variable, LAIT, GABA2, HDL, PLI, eMECA)

# Exportation en Excel ----
write.xlsx(results, file = "Table S1.xlsx")
