# Date : 13 novembre 2025
# Auteur : Jacob Lessard-Lord

# Projet PLC

# Extraction des résultats des statistiques des différentes étapes de validation

# Untargeted

# Set working directory ----
setwd("C:/Users/jackl/OneDrive/Bureau/École/Postdoctorat/Chercheur/Jean-Philippe/PLC/Processing_PLC/Untargeted/Final_Processing/Data")

# Chargement des packages requis ----
library(tidyverse)
library(openxlsx)

# Importation des données ----

# RCT - Steps 2-4
load("Metabo/Tableaux S7-S8-S9-S10 - Steps 2-3-4.rda")

tableau_lait_rct <- tableau_lait
tableau_fromage_rct <- tableau_fromage
tableau_beurre_rct <- tableau_beurre
tableau_dairy_fat_rct <- tableau_dairy_fat

rm(tableau_lait, tableau_fromage, tableau_beurre, tableau_dairy_fat)

# RCT - Dose-response - Step 5
load("Metabo/Tableaux S7-S8-S10 - Step 5.rda")

tableau_lait_dose <- tableau_lait
tableau_fromage_dose <- tableau_fromage
tableau_dairy_fat_dose <- tableau_dairy_fat

rm(tableau_lait, tableau_fromage, tableau_dairy_fat)

# Cohort - Step 6
load("Metabo/Tableaux S7-S9-S10 - Step 6.rda")

tableau_lait_cohort <- tableau_lait
tableau_beurre_cohort <- tableau_beurre
tableau_dairy_fat_cohort <- tableau_dairy_fat

rm(tableau_lait, tableau_beurre, tableau_dairy_fat)

# Merge des tableaux ----

# Lait 
tableau_lait <- full_join(x = full_join(x = tableau_lait_rct,
                                        y = tableau_lait_dose,
                                        by = "Candidate"),
                          y = tableau_lait_cohort,
                          by = "Candidate") %>%
  select(-ID)

rm(tableau_lait_rct, tableau_lait_dose, tableau_lait_cohort)

# Fromage
tableau_fromage <- full_join(x = tableau_fromage_rct,
                             y = tableau_fromage_dose,
                             by = "Candidate")

rm(tableau_fromage_rct, tableau_fromage_dose)

# Beurre
tableau_beurre <- full_join(x = tableau_beurre_rct,
                            y = tableau_beurre_cohort,
                            by = "Candidate") %>%
  select(-ID)

rm(tableau_beurre_rct, tableau_beurre_cohort)

# Dairy fat
tableau_dairy_fat <- full_join(x = full_join(x = tableau_dairy_fat_rct,
                                             y = tableau_dairy_fat_dose,
                                             by = "Candidate"),
                               y = tableau_dairy_fat_cohort,
                               by = "Candidate") %>%
  select(-ID)

rm(tableau_dairy_fat_rct, tableau_dairy_fat_dose, tableau_dairy_fat_cohort)


# Mise en forme des tableaux ----
tableau_lait <- tableau_lait %>% arrange(Step_6_p_value_ffq, Step5_p_value, Step4_p_adj, Step3_cheese_p_adj, Step2_p_adj)

tableau_fromage <- tableau_fromage %>% arrange(Step5_p_value, Step4_p_adj, Step3_milk_p_adj, Step2_p_adj)

tableau_beurre <- tableau_beurre %>% arrange(Step_6_p_value_ffq, Step3_cheese_p_adj, Step2_p_adj)

tableau_dairy_fat <- tableau_dairy_fat %>% arrange(Step_6_p_value_ffq, Step5_p_value, Step4_p_adj, Step3_cheese_p_adj)

# Remettre en forme les p-values correctement (pas de notation scientifique)

## pattern: single digit + e/E-04 (or -4)
pat <- "^(\\d)[eE]-0?4$"

## Fonction
fix_e04 <- function(x) {
  str_replace(as.character(x), pat, "0.000\\1")
}


## Changement dans les data.frames
tableau_lait <- tableau_lait %>%
  mutate(across(contains("_p_"), fix_e04))

tableau_fromage <- tableau_fromage %>%
  mutate(across(contains("_p_"), fix_e04))

tableau_beurre <- tableau_beurre %>%
  mutate(across(contains("_p_"), fix_e04))

tableau_dairy_fat <- tableau_dairy_fat %>%
  mutate(across(contains("_p_"), fix_e04))


# Exportation des résultats ----

# Créer un classeur Excel
wb <- createWorkbook()

# Ajouter une feuille pour chaque data.frame
addWorksheet(wb, "Milk")
writeData(wb, sheet = "Milk", x = tableau_lait)

addWorksheet(wb, "Cheese")
writeData(wb, sheet = "Cheese", x = tableau_fromage)

addWorksheet(wb, "Butter")
writeData(wb, sheet = "Butter", x = tableau_beurre)

addWorksheet(wb, "CFF_dairy")
writeData(wb, sheet = "CFF_dairy", x = tableau_dairy_fat)

# Enregistrer le classeur
saveWorkbook(wb, file = "Tables S7-S8-S11-S12.xlsx", overwrite = TRUE)
