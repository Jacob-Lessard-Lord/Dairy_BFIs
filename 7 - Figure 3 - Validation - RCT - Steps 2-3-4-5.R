# Date : 16 septembre 2025
# Auteur : Jacob Lessard-Lord

# Projet PLC

# Diagramme de Venn des différents biomarqueurs

# Untargeted

# Set working directory ----
setwd("C:/Users/jackl/OneDrive/Bureau/École/Postdoctorat/Chercheur/Jean-Philippe/PLC/Processing_PLC/Untargeted/Final_Processing/Data")

# Chargement des packages requis ----
library(tidyverse)
library(ggpubr)

# Importation des données ----
load("Metabo/Untargeted - ALL.rda")

# Importation et mise en forme des biomarqueurs ----
load("Metabo/Biomarqueurs post-dereplication.rda")

# Subset des données pour le lait, le fromage et le beurre ----

# Sélection des biomarqueurs potentiels seulement
data_ech <- data_ech[, c("Projet", "Sujet", "Projet_Sujet","Traitement", metabolite)]

# Normalisation
data_ech[,  metabolite] <- scale(log(data_ech[,  metabolite]))

# Lait #
data_lait <- subset(data_ech, subset = Projet %in% c("LAIT", "GABA2") & Traitement %in% c("Controle", "Lait") & Projet_Sujet != "LAIT_16") %>%
  mutate(Traitement = droplevels(Traitement))

# Fromage #
data_fromage <- subset(data_ech, subset = Projet %in% c("HDL", "GABA2") & Traitement %in% c("Controle", "Fromage")) %>%
  mutate(Traitement = droplevels(Traitement))

# Beurre #
data_beurre <- subset(data_ech, subset = Projet == "HDL" & Traitement %in% c("Controle", "Beurre")& !Projet_Sujet %in% c("HDL_24", "HDL_68")) %>%
  mutate(Traitement = droplevels(Traitement))

# PLI # 
data_pli <- subset(data_ech, subset = Projet == "PLI" & Traitement %in% c("Controle", "Produits laitiers"))

# Dose-response lait # 
data_dose_response_lait <- rbind(data_ech %>% filter(Projet == "HDL" & Traitement == "Controle") %>% mutate(Lait_groupe = "None\n(0mL)"),
                                 data_ech %>% filter(Projet == "PLI" & Traitement == "Produits laitiers") %>% mutate(Lait_groupe = "Moderate\n(375mL)"),
                                 data_ech %>% filter(Traitement == "Lait") %>% mutate(Lait_groupe = "High\n(541 to\n1315mL)")) %>%
  mutate(Lait_groupe = factor(Lait_groupe, levels = c("None\n(0mL)", "Moderate\n(375mL)", "High\n(541 to\n1315mL)")))

# Série de boxplots pour illustrer chacun des critères ----
graph_step2 <- ggboxplot(data = data_lait %>% select(Projet_Sujet, Traitement, RPLC_Pos_0.75_160.1302mz, RPLC_Pos_1.07_158.1170mz) %>%
                           pivot_longer(cols = c(RPLC_Pos_0.75_160.1302mz, RPLC_Pos_1.07_158.1170mz), names_to = "Biomarqueur", values_to = "Intensity") %>%
                           mutate(Biomarqueur = factor(case_match(Biomarqueur,
                                                                  "RPLC_Pos_0.75_160.1302mz" ~ "δ-Valerobetaine",
                                                                  "RPLC_Pos_1.07_158.1170mz" ~ "Homostachydrine"),
                                                       levels = c("δ-Valerobetaine", "Homostachydrine"))),
                         x = "Traitement", y = "Intensity", fill = "Traitement", palette = c("#505900", "#425ae1"),
                         facet.by = "Biomarqueur", nrow = 2, scales = "free_y",
                         panel.labs.font = list(size = 10.5, face = "bold", color = "white"), panel.labs.background = list(fill = "black"),
                         ylab = "Normalized abundance",
                         legend = "none", title = "Step 2\n(Plausibility)",
                         outlier.shape = NA) +
  geom_point(alpha = 0.25) +
  geom_line(aes(group = Projet_Sujet), alpha = 0.25) + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 11),
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 10),
        axis.title.x = element_blank()) +
  scale_x_discrete(labels = c("No\ndairy", "Milk"))

graph_step3_cheese <- ggboxplot(data = data_fromage %>% select(Projet_Sujet, Traitement, RPLC_Pos_0.75_160.1302mz, RPLC_Pos_1.07_158.1170mz) %>%
                                  pivot_longer(cols = c(RPLC_Pos_0.75_160.1302mz, RPLC_Pos_1.07_158.1170mz), names_to = "Biomarqueur", values_to = "Intensity") %>%
                                  mutate(Biomarqueur = factor(case_match(Biomarqueur,
                                                                         "RPLC_Pos_0.75_160.1302mz" ~ "δ-Valerobetaine",
                                                                         "RPLC_Pos_1.07_158.1170mz" ~ "Homostachydrine"),
                                                              levels = c("δ-Valerobetaine", "Homostachydrine"))),
                               x = "Traitement", y = "Intensity", fill = "Traitement", palette = c("#505900", "#fd533d"),
                               facet.by = "Biomarqueur", nrow = 2, scales = "free_y",
                               panel.labs.font = list(size = 10.5, face = "bold", color = "white"), panel.labs.background = list(fill = "black"),
                               ylab = "Normalized abundance",
                               legend = "none", title = "Step 3\n(Specificity)",
                               outlier.shape = NA) +
  geom_point(alpha = 0.25) +
  geom_line(aes(group = Projet_Sujet), alpha = 0.25) + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 11),
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 10),
        axis.title.x = element_blank()) +
  scale_x_discrete(labels = c("No\ndairy", "Cheese"))

graph_step3_butter <- ggboxplot(data = data_beurre %>% select(Projet_Sujet, Traitement, RPLC_Pos_0.75_160.1302mz, RPLC_Pos_1.07_158.1170mz) %>%
                                  pivot_longer(cols = c(RPLC_Pos_0.75_160.1302mz, RPLC_Pos_1.07_158.1170mz), names_to = "Biomarqueur", values_to = "Intensity") %>%
                                  mutate(Biomarqueur = factor(case_match(Biomarqueur,
                                                                         "RPLC_Pos_0.75_160.1302mz" ~ "δ-Valerobetaine",
                                                                         "RPLC_Pos_1.07_158.1170mz" ~ "Homostachydrine"),
                                                              levels = c("δ-Valerobetaine", "Homostachydrine"))),
                               x = "Traitement", y = "Intensity", fill = "Traitement", palette = c("#505900", "#FFD166"),
                               facet.by = "Biomarqueur", nrow = 2, scales = "free_y",
                               panel.labs.font = list(size = 10.5, face = "bold", color = "white"), panel.labs.background = list(fill = "black"),
                               ylab = "Normalized abundance",
                               legend = "none", title = "Step 3\n(Specificity)",
                               outlier.shape = NA) +
  geom_point(alpha = 0.25) +
  geom_line(aes(group = Projet_Sujet), alpha = 0.25) + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 11),
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 10),
        axis.title.x = element_blank()) +
  scale_x_discrete(labels = c("No\ndairy", "Butter"))

graph_step4 <- ggboxplot(data = data_pli %>% select(Projet_Sujet, Traitement, RPLC_Pos_0.75_160.1302mz, RPLC_Pos_1.07_158.1170mz) %>%
                           pivot_longer(cols = c(RPLC_Pos_0.75_160.1302mz, RPLC_Pos_1.07_158.1170mz), names_to = "Biomarqueur", values_to = "Intensity") %>%
                           mutate(Biomarqueur = factor(case_match(Biomarqueur,
                                                                  "RPLC_Pos_0.75_160.1302mz" ~ "δ-Valerobetaine",
                                                                  "RPLC_Pos_1.07_158.1170mz" ~ "Homostachydrine"),
                                                       levels = c("δ-Valerobetaine", "Homostachydrine"))),
                        x = "Traitement", y = "Intensity", fill = "Traitement", palette = c("#505900", "#02e0b1"),
                        facet.by = "Biomarqueur", nrow = 2, scales = "free_y",
                        panel.labs.font = list(size = 10.5, face = "bold", color = "white"), panel.labs.background = list(fill = "black"),
                        ylab = "Normalized abundance",
                        legend = "none", title = "Step 4\n(Robustness)",
                        outlier.shape = NA) +
  geom_point(alpha = 0.25) +
  geom_line(aes(group = Projet_Sujet), alpha = 0.25) + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 11),
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 10),
        axis.title.x = element_blank()) +
  scale_x_discrete(labels = c("No\ndairy", "Dairy\nproducts"))

graph_step5 <- ggboxplot(data = data_dose_response_lait %>% select(Projet_Sujet, Lait_groupe, RPLC_Pos_0.75_160.1302mz, RPLC_Pos_1.07_158.1170mz) %>%
                           pivot_longer(cols = c(RPLC_Pos_0.75_160.1302mz, RPLC_Pos_1.07_158.1170mz), names_to = "Biomarqueur", values_to = "Intensity") %>%
                           mutate(Biomarqueur = factor(case_match(Biomarqueur,
                                                                  "RPLC_Pos_0.75_160.1302mz" ~ "δ-Valerobetaine",
                                                                  "RPLC_Pos_1.07_158.1170mz" ~ "Homostachydrine"),
                                                       levels = c("δ-Valerobetaine", "Homostachydrine"))),
                        x = "Lait_groupe", y = "Intensity", fill = "Lait_groupe", palette = c("#01903d", "#507bff", "#8a335d"),
                        facet.by = "Biomarqueur", nrow = 2, scales = "free_y",
                        panel.labs.font = list(size = 10.5, face = "bold", color = "white"), panel.labs.background = list(fill = "black"),
                        ylab = "Scaled area", xlab = "Milk intake", title = "Step 5\n(Dose-response)",
                        outlier.shape = NA) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 11),
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 10),
        legend.position = "none") +
  geom_point(alpha = 0.25)

# Figure finale ----
setwd("C:/Users/jackl/OneDrive/Bureau/École/Postdoctorat/Chercheur/Jean-Philippe/PLC/Article/Biomarqueurs/Figures")

tiff(filename = "Figure 3 - Validation RCT.tiff", width = 25, height = 18, units = "cm", res = 600, compression = "lzw")

ggarrange(graph_step2, graph_step3_cheese, graph_step3_butter, graph_step4, graph_step5, ncol = 5, labels = c("A", "B", NA, "C", "D"),
          widths = c(1, 1, 1, 1, 1.25))

dev.off()

