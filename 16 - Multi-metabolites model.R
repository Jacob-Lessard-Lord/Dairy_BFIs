# Date : 11 septembre 2025
# Auteur : Jacob Lessard-Lord

# Projet PLC

# Entrainement du multi-metabolites model

# Set working directory ----
setwd("C:/Users/jackl/OneDrive/Bureau/École/Postdoctorat/Chercheur/Jean-Philippe/PLC/Processing_PLC/Untargeted/Final_Processing/Data")

# Chargement des packages requis ----
library(tidyverse)
library(ggpubr)
library(car)
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

# Sélection des biomarqueurs potentiels seulement
data_ech <- data_ech[, c("Projet", "Sujet", "Projet_Sujet","Traitement", all_biomarkers)]

# Scaling
data_ech[,  all_biomarkers] <- scale(log(data_ech[,  all_biomarkers]))

# Subset des données pour le lait et le fromage (pour entrainement du multi-metabolite model) ----

# Lait #
data_lait <- subset(data_ech, subset = Projet %in% c("LAIT", "GABA2") & Traitement %in% c("Controle", "Lait") & Projet_Sujet != "LAIT_16") %>%
  mutate(Traitement = droplevels(Traitement))

# Fromage #
data_fromage <- subset(data_ech, subset = Projet %in% c("HDL", "GABA2") & Traitement %in% c("Controle", "Fromage")) %>%
  mutate(Traitement = droplevels(Traitement))

# Beurre #
data_beurre <- subset(data_ech, subset = Projet == "HDL" & Traitement %in% c("Controle", "Beurre")) %>%
  mutate(Traitement = droplevels(Traitement))

# Dairy fat # 
data_dairy_fat <- subset(data_ech, subset = Projet == "HDL" & Traitement %in% c("Controle", "Fromage", "Beurre")) %>%
  mutate(Traitement = factor(case_match(Traitement,
                                        "Fromage" ~ "Dairy fat",
                                        "Beurre" ~ "Dairy fat",
                                        .default = Traitement)))

# Selection des données au baseline des RCT + eMECA ----
data_ech <- subset(data_ech, Traitement == "Baseline" | Projet == "eMECA")

# Modèles de régression logistique ----

# Lait # 
modele_lait <- glm(as.formula(paste("Traitement ~", paste(biomarqueur_lait, collapse = " + "))),
                   data = data_lait, family = "binomial")
summary(modele_lait)

# Verification de la multicolinéarite
vif(modele_lait) # VIF < 10, c'est OK


# Beurre # 
modele_beurre <- glm(as.formula(paste("Traitement ~", paste(biomarqueur_beurre, collapse = " + "))),
                   data = data_beurre, family = "binomial")
summary(modele_beurre)

# Verification de la multicolinéarite
vif(modele_beurre) # VIF < 10, c'est OK


# Dairy fat # 
modele_dairy_fat <- glm(as.formula(paste("Traitement ~", paste(biomarqueur_dairy_fat, collapse = " + "))),
                   data = data_dairy_fat, family = "binomial")
summary(modele_dairy_fat)

# Verification de la multicolinéarite
vif(modele_dairy_fat) # VIF < 10, c'est OK

# Prédiction du score métabolomique dans la cohorte ----
data_ech <- data_ech %>% 
  mutate(Score_lait = predict(modele_lait, newdata = ., type = "response"),
         Score_beurre = predict(modele_beurre, newdata = ., type = "response"),
         Score_dairy_fat = predict(modele_dairy_fat, newdata = ., type = "response"))

# Enregistrement des données
data_score <- data_ech %>% select(Projet, Sujet, Projet_Sujet, Score_lait, Score_beurre, Score_dairy_fat)

save(data_score, file = "Metabo/Score multi-metabolite.rda")
