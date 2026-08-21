# Date : 21 août 2026
# Auteur : Jacob Lessard-Lord

# Projet PLC

# ÉValuation de la performance des biomarqueurs en contexte observationnel

# Reproductibilité à court et long terme

# Caractéristiques des participants

# Set working directory ----
setwd("C:/Users/jackl/OneDrive/Bureau/École/Postdoctorat/Chercheur/Jean-Philippe/PLC/Processing_PLC/Untargeted/Final_Processing/Data")

# Chargement des packages requis ----
library(tidyverse)
library(rstatix)
library(openxlsx)

# Importation des données ----
load("FH/Data FFQ FH.rda")

# Mise en forme des données ----

# Sélection du baseline
data_ffq_fh_cac_lpa <- data_ffq_fh_cac_lpa %>% filter(Projet == "FH-LPA")
data_ffq_fh_diet <- data_ffq_fh_diet %>% filter(Visites == "V0") 


# Combinaison des datasets
data_repeated <- bind_rows(
  data_ffq_fh_cac_lpa %>% select(Projet, nomat, age, sex, bmi, Lait, Fromage, Beurre),
  data_ffq_fh_diet %>% select(Projet, nomat, age, sex, bmi, Lait, Fromage, Beurre)
)

# Calcul des statistiques descriptives ----
results <- rbind(data_repeated %>%
                   group_by(Projet) %>%
                   get_summary_stats(type = "mean_sd") %>%
                   mutate(Results_all = paste(round(mean, 1), " ± ", round(sd, 1))) %>%
                   select(-c(n, mean, sd)),
                 data_repeated %>% freq_table(Projet, sex) %>% mutate(Results_all = paste0(n, " (", prop, "%)")) %>% select(-c(n, prop)) %>% rename (variable = sex)) %>%
  pivot_wider(names_from = "Projet", values_from = "Results_all")

# Exportation en Excel ----
write.xlsx(results, file = "Table S14.xlsx")