# Date : 15 septembre 2025
# Auteur : Jacob Lessard-Lord

# Projet PLC

# Figure 5A

# Set working directory ----
setwd("C:/Users/jackl/OneDrive/Bureau/École/Postdoctorat/Chercheur/Jean-Philippe/PLC/Processing_PLC/Untargeted/Final_Processing/Data/Analyse aliments")

# Chargement des packages requis ----
library(tidyverse)
library(ggpubr)

# Importation des données ----
load("Food analysis.rda")

# Barplot ----
setwd("C:/Users/jackl/OneDrive/Bureau/École/Postdoctorat/Chercheur/Jean-Philippe/PLC/Article/Biomarqueurs/Figures")

tiff(filename = "Figure 5 - Barplot Aliment.tiff", width = 20, height = 25, units = "cm", res = 600, compression = "lzw")

ggbarplot(data = data_aliment %>% pivot_longer(cols = c("dVB", "Homostachydrine", "SM_35_1", "SM_43_1", "Cer_43_1", "SM_43_2", "PC_17_0_18_1"),
                                               names_to = "BFI", values_to = "Abundance") %>%
            mutate(BFI = factor(case_match(BFI,
                                           "dVB" ~ "δ-Valerobetaine",
                                           "Homostachydrine" ~ "Homostachydrine",
                                           "SM_35_1" ~ "SM d18:1/17:0",
                                           "SM_43_1" ~ "SM 43:1",
                                           "Cer_43_1" ~ "Cer d18:1/25:0",
                                           "SM_43_2" ~ "SM 43:2",
                                           "PC_17_0_18_1" ~ "PC (17:0/18:1)"),
                                levels = c("δ-Valerobetaine", "Homostachydrine", "SM d18:1/17:0", "SM 43:1", "Cer d18:1/25:0", "SM 43:2", "PC (17:0/18:1)"))),
          x = "Aliment", y = "Abundance", fill = "Food",
          palette = c("#425ae1", "#fd533d", "#06D6A0", "#ff93da", "#FFD166"),
          facet.by = "BFI", ncol = 4, nrow = 2,
          panel.labs.font = list(size = 10.5, face = "bold", color = "white"), panel.labs.background = list(fill = "black"),
          ylab = "Relative abundance (%)",
          legend = "bottom") +
  theme(axis.title.y = element_blank(),
        axis.title.x = element_text(size = 10),
        axis.text = element_text(size = 10),
        legend.text = element_text(size = 9),
        legend.title = element_text(size = 9),
        legend.box.spacing = unit(4, "pt"),
        legend.margin = margin(0, 0, 0, 0)) +
  coord_flip()

dev.off()
