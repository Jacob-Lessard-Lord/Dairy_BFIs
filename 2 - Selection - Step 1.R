# Date : 2 mai 2025
# Auteur : Jacob Lessard-Lord

# Projet PLC

# Identification de la signature metabolomique du lait et du fromage 

# Selection de variables par XGBoost

# Nouvelle approche d'Elsa

# Set seed (pour la reproductibilité) ----
set.seed(2608) 

# Set working directory ----
setwd("C:/Users/jackl/OneDrive/Bureau/École/Postdoctorat/Chercheur/Jean-Philippe/PLC/Processing_PLC/Untargeted/Final_Processing/Data")

# Chargement des packages requis ----
library(tidyverse)
library(caret)
library(xgboost)
library(ggpubr)
library(MLmetrics)
library(pROC)
library(rstatix)

# Importation des données ----
load("Metabo/Untargeted - ALL.rda")

# Subset des donnees pour le lait et le fromage ----
data_lait <- subset(data_ech,
                    subset = Projet %in% c("LAIT", "GABA2") & Traitement %in% c("Controle", "Lait") & Projet_Sujet != "LAIT_16",
                    select = c("Projet_Sujet", "Projet", "Traitement", metabolite)) %>%
  mutate(Traitement = droplevels(Traitement))

data_fromage <- subset(data_ech, subset = Projet %in% c("GABA2", "HDL") & Traitement %in% c("Controle", "Fromage"),
                       select = c("Projet_Sujet", "Projet", "Traitement", metabolite)) %>%
  mutate(Traitement = droplevels(Traitement))


data_beurre <- subset(data_ech, subset = Projet == "HDL" & Traitement %in% c("Controle", "Beurre") & !Projet_Sujet %in% c("HDL_24", "HDL_68"),
                       select = c("Projet_Sujet", "Projet", "Traitement", metabolite)) %>%
  mutate(Traitement = droplevels(Traitement))

data_lait_fromage <- subset(data_ech, subset = Projet == "GABA2" & Traitement %in% c("Lait", "Fromage"),
                      select = c("Projet_Sujet", "Projet", "Traitement", metabolite)) %>%
  mutate(Traitement = droplevels(Traitement))

data_fromage_beurre <- subset(data_ech, subset = Projet == "HDL" & Traitement %in% c("Fromage", "Beurre") & !Projet_Sujet %in% c("HDL_24", "HDL_68"),
                      select = c("Projet_Sujet", "Projet", "Traitement", metabolite)) %>%
  mutate(Traitement = droplevels(Traitement))

# Split des données pour cross-validation ----

## Nombre de splits désirés ##
nbr_split <- 25

# Lait
sujet_lait <- rbind(subset(data_ech, Projet %in% c("LAIT", "GABA2") & Traitement == "Lait", select = c("Projet_Sujet", "Projet", "Traitement")))
sujet_lait <- sujet_lait %>%
  mutate(Traitement = droplevels(Traitement))

## Partition pour les n splits
partition_lait <- createDataPartition(y = sujet_lait$Projet, p = 0.8, list = FALSE, times = nbr_split)
colnames(partition_lait) <- 1:nbr_split

# Fromage
sujet_fromage <- rbind(subset(data_ech, Projet %in% c("HDL", "GABA2") & Traitement == "Fromage", select = c("Projet_Sujet", "Projet", "Traitement")))
sujet_fromage <- sujet_fromage %>%
  mutate(Traitement = droplevels(Traitement))

# Partition pour les n splits
partition_fromage <- createDataPartition(y = sujet_fromage$Projet, p = 0.8, list = FALSE, times = nbr_split)
colnames(partition_fromage) <- 1:nbr_split

# Beurre
sujet_beurre <- rbind(subset(data_ech, Projet == "HDL" & Traitement == "Beurre" & !Projet_Sujet %in% c("HDL_24", "HDL_68"),
                             select = c("Projet_Sujet", "Projet", "Traitement")))
sujet_beurre <- sujet_beurre %>%
  mutate(Traitement = droplevels(Traitement))

# Partition pour les n splits
partition_beurre <- createDataPartition(y = sujet_beurre$Projet, p = 0.8, list = FALSE, times = nbr_split)
colnames(partition_beurre) <- 1:nbr_split

# Lait vs Fromage
sujet_lait_fromage <- rbind(subset(data_ech, Projet == "GABA2" & Traitement == "Fromage",
                                   select = c("Projet_Sujet", "Projet", "Traitement")))
sujet_lait_fromage <- sujet_lait_fromage %>%
  mutate(Traitement = droplevels(Traitement))

# Partition pour les n splits
partition_lait_fromage <- createDataPartition(y = sujet_lait_fromage$Projet, p = 0.8, list = FALSE, times = nbr_split)
colnames(partition_lait_fromage) <- 1:nbr_split

# Fromage vs Beurre
sujet_fromage_beurre <- rbind(subset(data_ech, Projet == "HDL" & Traitement == "Fromage" & !Projet_Sujet %in% c("HDL_24", "HDL_68"),
                                     select = c("Projet_Sujet", "Projet", "Traitement")))
sujet_fromage_beurre <- sujet_fromage_beurre %>%
  mutate(Traitement = droplevels(Traitement))

# Partition pour les n splits
partition_fromage_beurre <- createDataPartition(y = sujet_fromage_beurre$Projet, p = 0.8, list = FALSE, times = nbr_split)
colnames(partition_fromage_beurre) <- 1:nbr_split


# Controle pour machine learning ----
myControl_base <- trainControl(
  method = "cv", # Méthode de sélection par validation croisée
  number = 5, # Nombre d'échantillon de validation croisée
  summaryFunction = twoClassSummary, # Important pour obtenir des statistiques intéressantes pour les variables binaires. Retirer cette ligne pour les variables continues.
  #summaryFunction = fbeta,
  classProbs = T, # Important pour obtenir des statistiques intéressantes pour les variables binaires. Retirer cette ligne pour les variables continues
  verboseIter = FALSE # Permet de suivre l'avancement de l'algorithme.  Je vous suggère de le mettre à TRUE
)

# Fonction pour l'evaluation de modele ----
modele_evaluation <- function(modele, dataset, pos_class){
  
  set.seed(2608)
  
  # Prédiciton du modèle
  prediction <- predict(modele, newdata = as.matrix(subset(dataset, select = -Traitement)))
  prediction <- ifelse(test = prediction < 0.5, yes = "Controle", no = pos_class)
  
  # Matrice de confusion du modèle sur le dataset
  confusion <- MLmetrics::ConfusionMatrix( y_pred = prediction, y_true = dataset$Traitement)
  
  # Spécificité du modèle
  specificite <- MLmetrics::Specificity(y_pred = prediction, y_true = dataset$Traitement, positive = pos_class)
  
  # Sensibilite du modele
  sensibilite <- MLmetrics::Sensitivity(y_pred = prediction, y_true = dataset$Traitement, positive = pos_class)
  ### Analyse de la courbe ROC ###
  
  # Objet ROC
  roc_rf <- roc(response = factor(dataset$Traitement, ordered = T), predictor =  factor(prediction, ordered = TRUE))
  
  # Courbe ROC
  g_roc <- ggroc(roc_rf,  alpha = 0.5, colour = "red", linetype = 1, size = 2,
                 legacy.axes = TRUE) +
    theme_pubr() +
    ggtitle("Courbe ROC") +
    labs(x = "1 - Spécificité", y = "Sensibilité") +
    geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1),
                 linetype = 2) +
    theme(plot.title = element_text(hjust = 0.5))
  print(g_roc)
  
  # Aire sous la courbe ROC
  aire_courbe_roc <- pROC::auc(roc_rf)
  
  # Retour des objets
  valeurs <- list("matrice_confusion" = confusion, "valeur_specificite" = specificite,
                  "courbe_roc" =  g_roc, "auc" = aire_courbe_roc, "valeur_sensibilite"=sensibilite)
  
  return(valeurs)
}

# Validation croisée pour l'optimisation des hyperparamètres à chacun des n splits ----

## Lait ##

# Creation d'une liste avec les meilleurs hyperparametres par split
hyperparametre_lait <- list()

# Entrainement de n modeles
for (i in 1:nbr_split) {
  
  # Entrainement du modele
  modele <- train(
    form = Traitement ~ .,  
    data = subset(data_lait, subset = Projet_Sujet %in% sujet_lait[partition_lait[, i], ]$Projet_Sujet, select = c("Traitement", metabolite)),
    method = "xgbTree", # Le nom de la méthode
    trControl = myControl_base # Les contrôles créés avec trainControl
  )
  
  hyperparametre_lait[[paste0("Split_", i)]] <- modele$bestTune
  rm(modele)
}

## Fromage ##

# Creation d'une liste avec les meilleurs hyperparametres par split
hyperparametre_fromage <- list()

# Entrainement de n modeles
for (i in 1:nbr_split) {
  
  # Entrainement du modele
  modele <- train(
    form = Traitement ~ .,  
    data = subset(data_fromage, subset = Projet_Sujet %in% sujet_fromage[partition_fromage[, i], ]$Projet_Sujet, select = c("Traitement", metabolite)),
    method = "xgbTree", # Le nom de la méthode
    trControl = myControl_base # Les contrôles créés avec trainControl
  )
  
  hyperparametre_fromage[[paste0("Split_", i)]] <- modele$bestTune
  rm(modele)
}

## Beurre ##

# Creation d'une liste avec les meilleurs hyperparametres par split
hyperparametre_beurre <- list()

# Entrainement de n modeles
for (i in 1:nbr_split) {
  
  # Entrainement du modele
  modele <- train(
    form = Traitement ~ .,  
    data = subset(data_beurre, subset = Projet_Sujet %in% sujet_beurre[partition_beurre[, i], ]$Projet_Sujet, select = c("Traitement", metabolite)),
    method = "xgbTree", # Le nom de la méthode
    trControl = myControl_base # Les contrôles créés avec trainControl
  )
  
  hyperparametre_beurre[[paste0("Split_", i)]] <- modele$bestTune
  rm(modele)
}

## Lait vs Fromage ##

# Creation d'une liste avec les meilleurs hyperparametres par split
hyperparametre_lait_fromage <- list()

# Entrainement de n modeles
for (i in 1:nbr_split) {
  
  # Entrainement du modele
  modele <- train(
    form = Traitement ~ .,  
    data = subset(data_lait_fromage, subset = Projet_Sujet %in% sujet_lait_fromage[partition_lait_fromage[, i], ]$Projet_Sujet, select = c("Traitement", metabolite)),
    method = "xgbTree", # Le nom de la méthode
    trControl = myControl_base # Les contrôles créés avec trainControl
  )
  
  hyperparametre_lait_fromage[[paste0("Split_", i)]] <- modele$bestTune
  rm(modele)
}

## Fromage vs Beurre ##

# Creation d'une liste avec les meilleurs hyperparametres par split
hyperparametre_fromage_beurre <- list()

# Entrainement de n modeles
for (i in 1:nbr_split) {
  
  # Entrainement du modele
  modele <- train(
    form = Traitement ~ .,  
    data = subset(data_fromage_beurre, subset = Projet_Sujet %in% sujet_fromage_beurre[partition_fromage_beurre[, i], ]$Projet_Sujet, select = c("Traitement", metabolite)),
    method = "xgbTree", # Le nom de la méthode
    trControl = myControl_base # Les contrôles créés avec trainControl
  )
  
  hyperparametre_fromage_beurre[[paste0("Split_", i)]] <- modele$bestTune
  rm(modele)
}

# Sélection des features par processus itératif ----

# Retrait du feature le plus important après chaque itération de n splits de classification

## Lait ##
results_xgb_select_lait <- data.frame(AUC = NA,
                                      Specificite = NA,
                                      Sensibilite = NA,
                                      Split = NA)

metabolite_xgb_lait <- metabolite
met_xgb_lait <- vector()

set.seed(2608)
for (w in 1:1000) {
  # Initialisation des data.frames de resultats
  results_xgb <- data.frame(AUC = NA,
                            Specificite = NA,
                            Sensibilite = NA)
  
  # Initialisation data.frame importance des metabolites
  importance_metabolite <- data.frame(Metabolite = metabolite_xgb_lait)
  
  # Entrainement de n modeles
  for (i in 1:nbr_split) {
    
    # Entraintement du modele
    modele <- xgboost(data = as.matrix(subset(data_lait, subset = Projet_Sujet %in% sujet_lait[partition_lait[, i], ]$Projet_Sujet, select = metabolite_xgb_lait)),
                      label = ifelse(test = subset(data_lait, subset = Projet_Sujet %in% sujet_lait[partition_lait[, i], ]$Projet_Sujet)$Traitement == "Controle",
                                     yes = 0, no = 1),
                      objective = "binary:logistic",
                      nrounds = hyperparametre_lait[[paste0("Split_", i)]]$nrounds,
                      max_depth = hyperparametre_lait[[paste0("Split_", i)]]$max_depth,
                      eta = hyperparametre_lait[[paste0("Split_", i)]]$eta,
                      gamma = hyperparametre_lait[[paste0("Split_", i)]]$gamma,
                      colsamples_bytree = hyperparametre_lait[[paste0("Split_", i)]]$colsamples_bytree,
                      min_child_weight = hyperparametre_lait[[paste0("Split_", i)]]$min_child_weight,
                      subsample = hyperparametre_lait[[paste0("Split_", i)]]$subsample)
    
    # Validation du modele
    validation <- modele_evaluation(modele,
                                    subset(data_lait, subset = !Projet_Sujet %in% sujet_lait[partition_lait[, i], ]$Projet_Sujet,
                                           select = c("Traitement", metabolite_xgb_lait)),
                                    "Lait")
    
    # Extraction des metriques de validation
    results_xgb[i, "AUC"] <- validation$auc
    results_xgb[i, "Specificite"] <- validation$valeur_specificite
    results_xgb[i, "Sensibilite"] <- validation$valeur_sensibilite
    
    # Extraction des valeurs d'importance des metabolites
    xgb_importance <- xgb.importance(feature_names = metabolite_xgb_lait, model = modele)[, c("Feature", "Gain")]
    colnames(xgb_importance) <- c("Metabolite", "Importance")
    importance_metabolite <- merge(x = importance_metabolite, y = xgb_importance, by = "Metabolite", all = TRUE)
    colnames(importance_metabolite)[i + 1] <- paste0("Model_", i)
  }
  
  # Sauvegarde des resultats des metriques
  results_xgb$Split <- w
  results_xgb_select_lait <- rbind(results_xgb_select_lait, results_xgb)
  
  # Pondération de l'importance en fonction de l'AUC obtenu à chaque split
  for(i in 1:nbr_split) {
    importance_metabolite[, paste0("Model_", i)] <- ((importance_metabolite[, paste0("Model_", i)] * results_xgb[i, "AUC"]) / sum(results_xgb$AUC))
  }
  
  # Selection du meilleur feature
  importance_metabolite$Importance_somme <- apply(X = importance_metabolite %>% select(-Metabolite), MARGIN = 1, FUN = sum, na.rm = TRUE)
  met_xgb_lait <- c(met_xgb_lait, importance_metabolite[which.max(importance_metabolite$Importance_somme), "Metabolite"])
  
  # Retrait du metabolite le plus significatif
  metabolite_xgb_lait <- metabolite_xgb_lait[!metabolite_xgb_lait %in% importance_metabolite[which.max(importance_metabolite$Importance_somme), "Metabolite"]]
  
  # Arrêt de la boucle si limite inférieure de 95%CI de spécificité ou sensibilité est <= 0.5
  if(t.test(results_xgb$Sensibilite)$conf.int[1] <= 0.5 | t.test(results_xgb$Specificite)$conf.int[1] <= 0.5){
    break
  }
}

results_xgb_select_lait <- results_xgb_select_lait[-1, ]

## Fromage ##
results_xgb_select_fromage <- data.frame(AUC = NA,
                                         Specificite = NA,
                                         Sensibilite = NA,
                                         Split = NA)

metabolite_xgb_fromage <- metabolite
met_xgb_fromage <- vector()

set.seed(2608)
for (w in 1:1000) {
  # Initialisation des data.frames de resultats
  results_xgb <- data.frame(AUC = NA,
                            Specificite = NA,
                            Sensibilite = NA)
  
  # Initialisation data.frame importance des metabolites
  importance_metabolite <- data.frame(Metabolite = metabolite_xgb_fromage)
  
  # Entrainement de n modeles
  for (i in 1:nbr_split) {
    
    # Entraintement du modele
    modele <- xgboost(data = as.matrix(subset(data_fromage, subset = Projet_Sujet %in% sujet_fromage[partition_fromage[, i], ]$Projet_Sujet, select = metabolite_xgb_fromage)),
                      label = ifelse(test = subset(data_fromage, subset = Projet_Sujet %in% sujet_fromage[partition_fromage[, i], ]$Projet_Sujet)$Traitement == "Controle",
                                     yes = 0, no = 1),
                      objective = "binary:logistic",
                      nrounds = hyperparametre_fromage[[paste0("Split_", i)]]$nrounds,
                      max_depth = hyperparametre_fromage[[paste0("Split_", i)]]$max_depth,
                      eta = hyperparametre_fromage[[paste0("Split_", i)]]$eta,
                      gamma = hyperparametre_fromage[[paste0("Split_", i)]]$gamma,
                      colsamples_bytree = hyperparametre_fromage[[paste0("Split_", i)]]$colsamples_bytree,
                      min_child_weight = hyperparametre_fromage[[paste0("Split_", i)]]$min_child_weight,
                      subsample = hyperparametre_fromage[[paste0("Split_", i)]]$subsample)
    
    # Validation du modele
    validation <- modele_evaluation(modele,
                                    subset(data_fromage, subset = !Projet_Sujet %in% sujet_fromage[partition_fromage[, i], ]$Projet_Sujet,
                                           select = c("Traitement", metabolite_xgb_fromage)),
                                    "Fromage")
    
    # Extraction des metriques de validation
    results_xgb[i, "AUC"] <- validation$auc
    results_xgb[i, "Specificite"] <- validation$valeur_specificite
    results_xgb[i, "Sensibilite"] <- validation$valeur_sensibilite
    
    # Extraction des valeurs d'importance des metabolites
    xgb_importance <- xgb.importance(feature_names = metabolite_xgb_fromage, model = modele)[, c("Feature", "Gain")]
    colnames(xgb_importance) <- c("Metabolite", "Importance")
    importance_metabolite <- merge(x = importance_metabolite, y = xgb_importance, by = "Metabolite", all = TRUE)
    colnames(importance_metabolite)[i + 1] <- paste0("Model_", i)
  }
  
  # Sauvegarde des resultats des metriques
  results_xgb$Split <- w
  results_xgb_select_fromage <- rbind(results_xgb_select_fromage, results_xgb)
  
  # Pondération de l'importance en fonction de l'AUC obtenu à chaque split
  for(i in 1:nbr_split) {
    importance_metabolite[, paste0("Model_", i)] <- ((importance_metabolite[, paste0("Model_", i)] * results_xgb[i, "AUC"]) / sum(results_xgb$AUC))
  }
  
  # Selection du meilleur feature
  importance_metabolite$Importance_somme <- apply(X = importance_metabolite %>% select(-Metabolite), MARGIN = 1, FUN = sum, na.rm = TRUE)
  met_xgb_fromage <- c(met_xgb_fromage, importance_metabolite[which.max(importance_metabolite$Importance_somme), "Metabolite"])
  
  # Retrait du metabolite le plus significatif
  metabolite_xgb_fromage <- metabolite_xgb_fromage[!metabolite_xgb_fromage %in% importance_metabolite[which.max(importance_metabolite$Importance_somme), "Metabolite"]]
  
  # Arrêt de la boucle si limite inférieure de 95%CI de spécificité ou sensibilité est <= 0.5
  if(t.test(results_xgb$Sensibilite)$conf.int[1] <= 0.5 | t.test(results_xgb$Specificite)$conf.int[1] <= 0.5){
    break
  }
}

results_xgb_select_fromage <- results_xgb_select_fromage[-1, ]


## Beurre ##
results_xgb_select_beurre <- data.frame(AUC = NA,
                                         Specificite = NA,
                                         Sensibilite = NA,
                                         Split = NA)

metabolite_xgb_beurre <- metabolite
met_xgb_beurre <- vector()

set.seed(2608)
for (w in 1:1000) {
  # Initialisation des data.frames de resultats
  results_xgb <- data.frame(AUC = NA,
                            Specificite = NA,
                            Sensibilite = NA)
  
  # Initialisation data.frame importance des metabolites
  importance_metabolite <- data.frame(Metabolite = metabolite_xgb_beurre)
  
  # Entrainement de n modeles
  for (i in 1:nbr_split) {
    
    # Entraintement du modele
    modele <- xgboost(data = as.matrix(subset(data_beurre, subset = Projet_Sujet %in% sujet_beurre[partition_beurre[, i], ]$Projet_Sujet, select = metabolite_xgb_beurre)),
                      label = ifelse(test = subset(data_beurre, subset = Projet_Sujet %in% sujet_beurre[partition_beurre[, i], ]$Projet_Sujet)$Traitement == "Controle",
                                     yes = 0, no = 1),
                      objective = "binary:logistic",
                      nrounds = hyperparametre_beurre[[paste0("Split_", i)]]$nrounds,
                      max_depth = hyperparametre_beurre[[paste0("Split_", i)]]$max_depth,
                      eta = hyperparametre_beurre[[paste0("Split_", i)]]$eta,
                      gamma = hyperparametre_beurre[[paste0("Split_", i)]]$gamma,
                      colsamples_bytree = hyperparametre_beurre[[paste0("Split_", i)]]$colsamples_bytree,
                      min_child_weight = hyperparametre_beurre[[paste0("Split_", i)]]$min_child_weight,
                      subsample = hyperparametre_beurre[[paste0("Split_", i)]]$subsample)
    
    # Validation du modele
    validation <- modele_evaluation(modele,
                                    subset(data_beurre, subset = !Projet_Sujet %in% sujet_beurre[partition_beurre[, i], ]$Projet_Sujet,
                                           select = c("Traitement", metabolite_xgb_beurre)),
                                    "Beurre")
    
    # Extraction des metriques de validation
    results_xgb[i, "AUC"] <- validation$auc
    results_xgb[i, "Specificite"] <- validation$valeur_specificite
    results_xgb[i, "Sensibilite"] <- validation$valeur_sensibilite
    
    # Extraction des valeurs d'importance des metabolites
    xgb_importance <- xgb.importance(feature_names = metabolite_xgb_beurre, model = modele)[, c("Feature", "Gain")]
    colnames(xgb_importance) <- c("Metabolite", "Importance")
    importance_metabolite <- merge(x = importance_metabolite, y = xgb_importance, by = "Metabolite", all = TRUE)
    colnames(importance_metabolite)[i + 1] <- paste0("Model_", i)
  }
  
  # Sauvegarde des resultats des metriques
  results_xgb$Split <- w
  results_xgb_select_beurre <- rbind(results_xgb_select_beurre, results_xgb)
  
  # Pondération de l'importance en fonction de l'AUC obtenu à chaque split
  for(i in 1:nbr_split) {
    importance_metabolite[, paste0("Model_", i)] <- ((importance_metabolite[, paste0("Model_", i)] * results_xgb[i, "AUC"]) / sum(results_xgb$AUC))
  }
  
  # Selection du meilleur feature
  importance_metabolite$Importance_somme <- apply(X = importance_metabolite %>% select(-Metabolite), MARGIN = 1, FUN = sum, na.rm = TRUE)
  met_xgb_beurre <- c(met_xgb_beurre, importance_metabolite[which.max(importance_metabolite$Importance_somme), "Metabolite"])
  
  # Retrait du metabolite le plus significatif
  metabolite_xgb_beurre <- metabolite_xgb_beurre[!metabolite_xgb_beurre %in% importance_metabolite[which.max(importance_metabolite$Importance_somme), "Metabolite"]]
  
  # Arrêt de la boucle si limite inférieure de 95%CI de spécificité ou sensibilité est <= 0.5
  if(t.test(results_xgb$Sensibilite)$conf.int[1] <= 0.5 | t.test(results_xgb$Specificite)$conf.int[1] <= 0.5){
    break
  }
}

results_xgb_select_beurre <- results_xgb_select_beurre[-1, ]

## Lait vs Fromage ##
results_xgb_select_lait_fromage <- data.frame(AUC = NA,
                                                Specificite = NA,
                                                Sensibilite = NA,
                                                Split = NA)

metabolite_xgb_lait_fromage <- metabolite
met_xgb_lait_fromage <- vector()

set.seed(2608)
for (w in 1:1000) {
  # Initialisation des data.frames de resultats
  results_xgb <- data.frame(AUC = NA,
                            Specificite = NA,
                            Sensibilite = NA)
  
  # Initialisation data.frame importance des metabolites
  importance_metabolite <- data.frame(Metabolite = metabolite_xgb_lait_fromage)
  
  # Entrainement de n modeles
  for (i in 1:nbr_split) {
    
    # Entraintement du modele
    modele <- xgboost(data = as.matrix(subset(data_lait_fromage, subset = Projet_Sujet %in% sujet_lait_fromage[partition_lait_fromage[, i], ]$Projet_Sujet, select = metabolite_xgb_lait_fromage)),
                      label = ifelse(test = subset(data_lait_fromage, subset = Projet_Sujet %in% sujet_lait_fromage[partition_lait_fromage[, i], ]$Projet_Sujet)$Traitement == "Fromage",
                                     yes = 0, no = 1),
                      objective = "binary:logistic",
                      nrounds = hyperparametre_lait_fromage[[paste0("Split_", i)]]$nrounds,
                      max_depth = hyperparametre_lait_fromage[[paste0("Split_", i)]]$max_depth,
                      eta = hyperparametre_lait_fromage[[paste0("Split_", i)]]$eta,
                      gamma = hyperparametre_lait_fromage[[paste0("Split_", i)]]$gamma,
                      colsamples_bytree = hyperparametre_lait_fromage[[paste0("Split_", i)]]$colsamples_bytree,
                      min_child_weight = hyperparametre_lait_fromage[[paste0("Split_", i)]]$min_child_weight,
                      subsample = hyperparametre_lait_fromage[[paste0("Split_", i)]]$subsample)
    
    # Validation du modele
    validation <- modele_evaluation(modele,
                                    subset(data_lait_fromage, subset = !Projet_Sujet %in% sujet_lait_fromage[partition_lait_fromage[, i], ]$Projet_Sujet,
                                           select = c("Traitement", metabolite_xgb_lait_fromage)),
                                    "Lait")
    
    # Extraction des metriques de validation
    results_xgb[i, "AUC"] <- validation$auc
    results_xgb[i, "Specificite"] <- validation$valeur_specificite
    results_xgb[i, "Sensibilite"] <- validation$valeur_sensibilite
    
    # Extraction des valeurs d'importance des metabolites
    xgb_importance <- xgb.importance(feature_names = metabolite_xgb_lait_fromage, model = modele)[, c("Feature", "Gain")]
    colnames(xgb_importance) <- c("Metabolite", "Importance")
    importance_metabolite <- merge(x = importance_metabolite, y = xgb_importance, by = "Metabolite", all = TRUE)
    colnames(importance_metabolite)[i + 1] <- paste0("Model_", i)
  }
  
  # Sauvegarde des resultats des metriques
  results_xgb$Split <- w
  results_xgb_select_lait_fromage <- rbind(results_xgb_select_lait_fromage, results_xgb)
  
  # Pondération de l'importance en fonction de l'AUC obtenu à chaque split
  for(i in 1:nbr_split) {
    importance_metabolite[, paste0("Model_", i)] <- ((importance_metabolite[, paste0("Model_", i)] * results_xgb[i, "AUC"]) / sum(results_xgb$AUC))
  }
  
  # Selection du meilleur feature
  importance_metabolite$Importance_somme <- apply(X = importance_metabolite %>% select(-Metabolite), MARGIN = 1, FUN = sum, na.rm = TRUE)
  met_xgb_lait_fromage <- c(met_xgb_lait_fromage, importance_metabolite[which.max(importance_metabolite$Importance_somme), "Metabolite"])
  
  # Retrait du metabolite le plus significatif
  metabolite_xgb_lait_fromage <- metabolite_xgb_lait_fromage[!metabolite_xgb_lait_fromage %in% importance_metabolite[which.max(importance_metabolite$Importance_somme), "Metabolite"]]
  
  # Arrêt de la boucle si limite inférieure de 95%CI de spécificité ou sensibilité est <= 0.5
  if(t.test(results_xgb$Sensibilite)$conf.int[1] <= 0.5 | t.test(results_xgb$Specificite)$conf.int[1] <= 0.5){
    break
  }
}

results_xgb_select_lait_fromage <- results_xgb_select_lait_fromage[-1, ]


## Fromage vs Beurre ##
results_xgb_select_fromage_beurre <- data.frame(AUC = NA,
                                        Specificite = NA,
                                        Sensibilite = NA,
                                        Split = NA)

metabolite_xgb_fromage_beurre <- metabolite
met_xgb_fromage_beurre <- vector()

set.seed(2608)
for (w in 1:1000) {
  # Initialisation des data.frames de resultats
  results_xgb <- data.frame(AUC = NA,
                            Specificite = NA,
                            Sensibilite = NA)
  
  # Initialisation data.frame importance des metabolites
  importance_metabolite <- data.frame(Metabolite = metabolite_xgb_fromage_beurre)
  
  # Entrainement de n modeles
  for (i in 1:nbr_split) {
    
    # Entraintement du modele
    modele <- xgboost(data = as.matrix(subset(data_fromage_beurre, subset = Projet_Sujet %in% sujet_fromage_beurre[partition_fromage_beurre[, i], ]$Projet_Sujet, select = metabolite_xgb_fromage_beurre)),
                      label = ifelse(test = subset(data_fromage_beurre, subset = Projet_Sujet %in% sujet_fromage_beurre[partition_fromage_beurre[, i], ]$Projet_Sujet)$Traitement == "Fromage",
                                     yes = 0, no = 1),
                      objective = "binary:logistic",
                      nrounds = hyperparametre_fromage_beurre[[paste0("Split_", i)]]$nrounds,
                      max_depth = hyperparametre_fromage_beurre[[paste0("Split_", i)]]$max_depth,
                      eta = hyperparametre_fromage_beurre[[paste0("Split_", i)]]$eta,
                      gamma = hyperparametre_fromage_beurre[[paste0("Split_", i)]]$gamma,
                      colsamples_bytree = hyperparametre_fromage_beurre[[paste0("Split_", i)]]$colsamples_bytree,
                      min_child_weight = hyperparametre_fromage_beurre[[paste0("Split_", i)]]$min_child_weight,
                      subsample = hyperparametre_fromage_beurre[[paste0("Split_", i)]]$subsample)
    
    # Validation du modele
    validation <- modele_evaluation(modele,
                                    subset(data_fromage_beurre, subset = !Projet_Sujet %in% sujet_fromage_beurre[partition_fromage_beurre[, i], ]$Projet_Sujet,
                                           select = c("Traitement", metabolite_xgb_fromage_beurre)),
                                    "Beurre")
    
    # Extraction des metriques de validation
    results_xgb[i, "AUC"] <- validation$auc
    results_xgb[i, "Specificite"] <- validation$valeur_specificite
    results_xgb[i, "Sensibilite"] <- validation$valeur_sensibilite
    
    # Extraction des valeurs d'importance des metabolites
    xgb_importance <- xgb.importance(feature_names = metabolite_xgb_fromage_beurre, model = modele)[, c("Feature", "Gain")]
    colnames(xgb_importance) <- c("Metabolite", "Importance")
    importance_metabolite <- merge(x = importance_metabolite, y = xgb_importance, by = "Metabolite", all = TRUE)
    colnames(importance_metabolite)[i + 1] <- paste0("Model_", i)
  }
  
  # Sauvegarde des resultats des metriques
  results_xgb$Split <- w
  results_xgb_select_fromage_beurre <- rbind(results_xgb_select_fromage_beurre, results_xgb)
  
  # Pondération de l'importance en fonction de l'AUC obtenu à chaque split
  for(i in 1:nbr_split) {
    importance_metabolite[, paste0("Model_", i)] <- ((importance_metabolite[, paste0("Model_", i)] * results_xgb[i, "AUC"]) / sum(results_xgb$AUC))
  }
  
  # Selection du meilleur feature
  importance_metabolite$Importance_somme <- apply(X = importance_metabolite %>% select(-Metabolite), MARGIN = 1, FUN = sum, na.rm = TRUE)
  met_xgb_fromage_beurre <- c(met_xgb_fromage_beurre, importance_metabolite[which.max(importance_metabolite$Importance_somme), "Metabolite"])
  
  # Retrait du metabolite le plus significatif
  metabolite_xgb_fromage_beurre <- metabolite_xgb_fromage_beurre[!metabolite_xgb_fromage_beurre %in% importance_metabolite[which.max(importance_metabolite$Importance_somme), "Metabolite"]]
  
  # Arrêt de la boucle si limite inférieure de 95%CI de spécificité ou sensibilité est <= 0.5
  if(t.test(results_xgb$Sensibilite)$conf.int[1] <= 0.5 | t.test(results_xgb$Specificite)$conf.int[1] <= 0.5){
    break
  }
}

results_xgb_select_fromage_beurre <- results_xgb_select_fromage_beurre[-1, ]

# Sauvegarde des résultats ----
metabolite_lait <- met_xgb_lait
metabolite_fromage <- met_xgb_fromage
metabolite_beurre <- met_xgb_beurre

metabolite_lait_fromage <- met_xgb_lait_fromage
metabolite_fromage_beurre <- met_xgb_fromage_beurre

save(results_xgb_select_lait, metabolite_lait, results_xgb_select_fromage, metabolite_fromage, results_xgb_select_beurre, metabolite_beurre,
     results_xgb_select_lait_fromage, metabolite_lait_fromage, results_xgb_select_fromage_beurre, metabolite_fromage_beurre,
     file = "XGB - Selection variables.rda")
