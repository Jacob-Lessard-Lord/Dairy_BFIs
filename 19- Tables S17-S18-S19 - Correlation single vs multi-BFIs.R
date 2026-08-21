# Date : 21 août 2026
# Auteur : Jacob Lessard-Lord

# Projet PLC

# Corrélation entre chacun des biomarqueurs et le modèle multi-BFIs 

# Set working directory ----
setwd("C:/Users/jackl/OneDrive/Bureau/École/Postdoctorat/Chercheur/Jean-Philippe/PLC/Processing_PLC/Untargeted/Final_Processing/Data")

# Chargement des packages requis ----
library(tidyverse)
library(psych)
library(readxl)
library(writexl)

# Importation des données métabolomiques ----
load("Metabo/Untargeted - ALL.rda")

# Importation  des biomarqueurs ----
load("Metabo/Biomarqueurs finaux.rda")

# Importation des scores multi-métabolites ----
load("Metabo/Score multi-metabolite.rda")

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

# Ajout des scores
identification <- rbind(identification,
                        data.frame(Identification = c("Milk multi-BFIs model", "Butter multi-BFIs model", "CFF dairy multi-BFIs model"),
                                   Biomarqueur = c("Score_lait", "Score_beurre", "Score_dairy_fat")))
# Mise en forme des données ----

# Création d'un vecteur avec tous les biomarqueurs
all_biomarkers <- c(biomarqueur_lait, biomarqueur_beurre, biomarqueur_dairy_fat)

# Selection des données au baseline des RCT + eMECA
data_ech <- subset(data_ech, Traitement == "Baseline" | Projet == "eMECA", select = c("Projet", "Sujet", "Projet_Sujet", all_biomarkers))

# Scaling
data_ech[,  all_biomarkers] <- scale(log(data_ech[,  all_biomarkers]))

# Merge avec les données de score multi-metabolites
data_ech <- merge(x = data_score, data_ech, by = c("Projet", "Sujet", "Projet_Sujet"))


# Correlation de Spearman ----

# Lait #

# Select variables
data_lait <- data_ech %>%
  rename(!!!setNames(identification$Biomarqueur, identification$Identification)) %>%
  select(all_of(c("δ-Valerobetaine", "Homostachydrine", "Milk multi-BFIs model")))

# Spearman correlations
corr_lait <- corr.test(data_lait,
                       method = "spearman",
                       use = "pairwise",
                       adjust = "none")

# Create formatted table
corr_lait_table <- matrix(
  paste0("r = ", sprintf("%.2f", corr_lait$r),
    "; p = ", sprintf("%.4f", corr_lait$p)),
  nrow = nrow(corr_lait$r),
  dimnames = dimnames(corr_lait$r))

# Beurre #

# Select variables
data_beurre <- data_ech %>%
  rename(!!!setNames(identification$Biomarqueur, identification$Identification)) %>%
  select(all_of(c("PE (P-17:0/20:4, O-17:1/20:4)", "PC (P-17:0/20:4, O-17:1/20:4)", "Butter multi-BFIs model")))

# Spearman correlations
corr_beurre <- corr.test(data_beurre,
                       method = "spearman",
                       use = "pairwise",
                       adjust = "none")

# Create formatted table
corr_beurre_table <- matrix(
  paste0("r = ", sprintf("%.2f", corr_beurre$r),
         "; p = ", sprintf("%.4f", corr_beurre$p)),
  nrow = nrow(corr_beurre$r),
  dimnames = dimnames(corr_beurre$r))

# CFF dairy #

# Select variables
data_cff <- data_ech %>%
  rename(!!!setNames(identification$Biomarqueur, identification$Identification)) %>%
  select(all_of(c("SM d18:1/17:0", "PC 17:0/18:1", "SM 43:2", "SM 43:1", "Cer d18:1/25:0", "CFF dairy multi-BFIs model")))

# Spearman correlations
corr_cff <- corr.test(data_cff,
                       method = "spearman",
                       use = "pairwise",
                       adjust = "none")

# Create formatted table
corr_cff_table <- matrix(
  paste0("r = ", sprintf("%.2f", corr_cff$r),
         "; p = ", sprintf("%.4f", corr_cff$p)),
  nrow = nrow(corr_cff$r),
  dimnames = dimnames(corr_cff$r))



# Regrouper les tableaux dans une liste
correlation_matrices <- list(
  Milk = corr_lait_table %>%
    as.data.frame() %>%
    rownames_to_column("BFIs"),
  
  Butter = corr_beurre_table %>%
    as.data.frame() %>%
    rownames_to_column("BFIs"),
  
  CFF_dairy = corr_cff_table %>%
    as.data.frame() %>%
    rownames_to_column("BFIs")
)

correlation_matrices <- list(Milk = corr_lait_table,
                             Butter = corr_beurre_table,
                             CFF_dairy = corr_cff_table)

write_xlsx(correlation_matrices, path = "Table S17-S18-S19.xlsx")
