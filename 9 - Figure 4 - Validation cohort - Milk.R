# Date : 24 juillet 2025
# Auteur : Jacob Lessard-Lord

# Projet PLC

# ÉValuation de la performance des biomarqueurs en contexte observationnel

# Figure 4

# Set working directory ----
setwd("C:/Users/jackl/OneDrive/Bureau/École/Postdoctorat/Chercheur/Jean-Philippe/PLC/Processing_PLC/Untargeted/Final_Processing/Data")

# Chargement des packages requis ----
library(tidyverse)
library(ggpubr)
library(readxl)

# Importation des données ----
load("Metabo/Results Step 6.rda")

# Importation des identifications ----

# Lait
id_lait <- read_excel("Metabo/ID Biomarqueurs.xlsx", sheet = "Lait") %>% select(c(Chromatography, Polarity, RT, mz, Identification)) %>%
  mutate(ID = paste(Chromatography, Polarity, RT, mz, sep = "_"),
         label = if_else(is.na(Identification) | Identification %in% c("NA",""),
                         ID, Identification),
         # optional: wrap long labels so they fit
         label = str_wrap(label, width = 20)) %>%
  transmute(ID, label) %>%
  tibble::deframe()

# Graphique ----

# Lait #
legend_map <- c(
  "Q1 - Very low consumption (<15g)" = "<15 g/day",
  "Q2 - Low consumption (15-250g)"  = "15–250 g/day",
  "Q3 - High consumption (>250 g)"   = ">250 g/day"
)

graph_lait <- ggarrange(
  lait_ffq$plot + ggtitle("FFQ\n(n=342)") + theme(plot.title = element_text(hjust = 0.5, face = "bold"),
                                                  legend.text = element_text(size = 12),
                                                  legend.spacing.x  = unit(2, "pt"),
                                                  legend.spacing.y  = unit(0, "pt"),
                                                  legend.key.width  = unit(8, "pt"),
                                                  legend.key.height = unit(8, "pt"),
                                                  axis.title.x = element_blank()) + guides(color = guide_legend(ncol = 3)) + scale_x_discrete(labels = id_lait) +
    scale_color_manual(breaks = names(legend_map), labels = unname(legend_map), values = c("#01903d", "#507bff", "#8a335d")),
  
  lait_r24w_mean$plot + ggtitle("Mean of\n24h recalls\n(n=192)") + theme(plot.title = element_text(hjust = 0.5, face = "bold"),
                                                                         axis.title.x = element_blank()) + scale_x_discrete(labels = id_lait) +
    scale_color_manual(breaks = names(legend_map), labels = unname(legend_map), values = c("#01903d", "#507bff", "#8a335d")),
  
  lait_r24w_q3$plot + ggtitle("Day-of-draw\n24h recall\n(n=178)") + theme(plot.title = element_text(hjust = 0.5, face = "bold"),
                                                                          axis.title.x = element_blank()) + scale_x_discrete(labels = id_lait) +
    scale_color_manual(breaks = names(legend_map), labels = unname(legend_map), values = c("#01903d", "#507bff", "#8a335d")),
  
  ncol = 3, common.legend = TRUE, legend = "bottom", labels = "AUTO", widths = c(1.5, 1, 1)
)

# Exportation de la figure ----
setwd("C:/Users/jackl/OneDrive/Bureau/École/Postdoctorat/Chercheur/Jean-Philippe/PLC/Article/Biomarqueurs/Figures")

tiff(filename = "Figure 4 - Validation - Cohort - Milk.tiff", width = 21.5, height = 12.5, units = "cm", res = 600, compression = "lzw")
graph_lait
dev.off()