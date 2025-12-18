# Date : 24 juillet 2025
# Auteur : Jacob Lessard-Lord

# Projet PLC

# ÉValuation de la performance des biomarqueurs en contexte observationnel

# Figure S2

# Set working directory ----
setwd("C:/Users/jackl/OneDrive/Bureau/École/Postdoctorat/Chercheur/Jean-Philippe/PLC/Processing_PLC/Untargeted/Final_Processing/Data")

# Chargement des packages requis ----
library(tidyverse)
library(ggpubr)
library(readxl)

# Importation des données ----
load("Metabo/Results Step 6.rda")

# Importation des identifications ----

# Beurre
id_beurre <- read_excel("Metabo/ID Biomarqueurs.xlsx", sheet = "Beurre") %>% select(c(Chromatography, Polarity, RT, mz, Identification)) %>%
  mutate(ID = paste(Chromatography, Polarity, RT, mz, sep = "_"),
         label = if_else(is.na(Identification) | Identification %in% c("NA",""),
                         ID, Identification),
         # optional: wrap long labels so they fit
         label = str_wrap(label, width = 10)) %>%
  transmute(ID, label) %>%
  tibble::deframe()

# Dairy fat
id_dairy_fat <- read_excel("Metabo/ID Biomarqueurs.xlsx", sheet = "Dairy_fat") %>% select(c(Chromatography, Polarity, RT, mz, Identification)) %>%
  mutate(ID = paste(Chromatography, Polarity, RT, mz, sep = "_"),
         label = if_else(is.na(Identification) | Identification %in% c("NA",""),
                         ID, Identification),
         # optional: wrap long labels so they fit
         label = str_wrap(label, width = 10)) %>%
  transmute(ID, label) %>%
  tibble::deframe()

# Graphique ----
# Beurre #
legend_map <- c(
  "Q1 - No consumption (0g)" = "0 g/day",
  "Q2 - Low consumption (0-5g)"  = "0–5 g/day",
  "Q3 - High consumption (>5 g)"   = ">5 g/day"
)

graph_beurre <- ggarrange(
  beurre_ffq$plot + ggtitle("FFQ\n(n=342)") +  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
                                                     legend.text = element_text(size = 12),
                                                     legend.spacing.x  = unit(2, "pt"),
                                                     legend.spacing.y  = unit(0, "pt"),
                                                     legend.key.width  = unit(8, "pt"),
                                                     legend.key.height = unit(8, "pt"),
                                                     axis.title.x = element_blank()) + guides(color = guide_legend(ncol = 3)) + scale_x_discrete(labels = id_beurre) +
    scale_color_manual(breaks = names(legend_map), labels = unname(legend_map), values = c("#01903d", "#507bff", "#8a335d")),
  
  beurre_r24w_mean$plot + ggtitle("Mean of 24h recalls\n(n=192)") + theme(plot.title = element_text(hjust = 0.5, face = "bold"),
                                                                          axis.title.x = element_blank()) + scale_x_discrete(labels = id_beurre) +
    scale_color_manual(breaks = names(legend_map), labels = unname(legend_map), values = c("#01903d", "#507bff", "#8a335d")),
  
  beurre_r24w_q3$plot + ggtitle("Day-of-draw 24h recall\n(n=178)") + theme(plot.title = element_text(hjust = 0.5, face = "bold"),
                                                                           axis.title.x = element_blank()) + scale_x_discrete(labels = id_beurre) +
    scale_color_manual(breaks = names(legend_map), labels = unname(legend_map), values = c("#01903d", "#507bff", "#8a335d")),
  
  ncol = 3, common.legend = TRUE, legend = "bottom"
)

# Dairy fat # 
legend_map <- c(
  "Q1 - Very low consumption (<15g)" = "<15 g/day",
  "Q2 - Low consumption (15-60g)"  = "15-60 g/day",
  "Q3 - High consumption (>60 g)"   = ">60 g/day"
)

graph_dairy_fat <- ggarrange(
  dairy_fat_ffq$plot + theme(plot.title = element_text(hjust = 0.5, face = "bold"),
                             legend.text = element_text(size = 12),
                             legend.spacing.x  = unit(2, "pt"),
                             legend.spacing.y  = unit(0, "pt"),
                             legend.key.width  = unit(8, "pt"),
                             legend.key.height = unit(8, "pt"),
                             axis.title.x = element_blank(),
                             legend.position = "right") + guides(color = guide_legend(ncol = 3)) + scale_x_discrete(labels = id_dairy_fat) +
    scale_color_manual(breaks = names(legend_map), labels = unname(legend_map), values = c("#01903d", "#507bff", "#8a335d")),
  
  dairy_fat_r24w_mean$plot + theme(plot.title = element_text(hjust = 0.5, face = "bold"),
                                   axis.title.x = element_blank()) + scale_x_discrete(labels = id_dairy_fat) +
    scale_color_manual(breaks = names(legend_map), labels = unname(legend_map), values = c("#01903d", "#507bff", "#8a335d")),
  
  dairy_fat_r24w_q3$plot + theme(plot.title = element_text(hjust = 0.5, face = "bold"),
                                 axis.title.x = element_blank()) + scale_x_discrete(labels = id_dairy_fat) +
    scale_color_manual(breaks = names(legend_map), labels = unname(legend_map), values = c("#01903d", "#507bff", "#8a335d")),
  
  ncol = 3, common.legend = TRUE, legend = "bottom"
)

# Annotation
add_top_label <- function(p, lab) {
  annotate_figure(
    p,
    top = text_grob(lab, x = 0, hjust = 0, face = "bold", size = 14)
  )
}

graph_dairy_fat <- add_top_label(graph_dairy_fat,   "B) CFF dairy")

# Exportation de la figure ----
setwd("C:/Users/jackl/OneDrive/Bureau/École/Postdoctorat/Chercheur/Jean-Philippe/PLC/Article/Biomarqueurs/Figures")

tiff(filename = "Figure S2 - Validation - Cohort - Other.tiff", width = 32.5*1.15, height = 20*1.15, units = "cm", res = 600, compression = "lzw")
ggarrange(graph_beurre, graph_dairy_fat, nrow = 2,
          labels = c("A) Butter", NA),
          label.x = 0, label.y = 1,   # top-left anchor
          hjust = 0, vjust = 2       # left-align multi-line labels
)
dev.off()
