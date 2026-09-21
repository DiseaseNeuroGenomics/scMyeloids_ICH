# ============================================================================
# Figure1_C_umap.R
# Project: Decoding Immune Landscapes and Therapeutic Targets in Intracerebral Hemorrhage
# Author:  Dimitrios Kyriakis
# ============================================================================

# ==============================================================================
# Figure 1, Panel C — Immune cell UMAP colored by CellType
#
# Input:  2024_11_Final_Immune.rds
# Output: FIG_DIR/Figure1_C_umap.(pdf|png)
# ==============================================================================

source("../../00_setup.R")
source("../../00_helpers_stats.R")
source("../../00_theme_colors.R")

Immune <- readRDS(file.path(BASE_DIR, "2024_11_Final_Immune.rds"))

p_C_umap <- DimPlot(
  Immune,
  pt.size    = 0.01,
  cols       = Project_Colors,
  group.by   = "CellType",
  reduction  = "harmony_umap",
  raster     = FALSE
) +
  ggtitle("") +
  theme(
    plot.title = element_text(size = 12, face = "bold"),
    legend.title = element_text(size = 10),
    legend.text  = element_text(size = 11),
    legend.position = c(0.98, 0.38),
    legend.justification = c("right", "top"),
    legend.background = element_rect(fill = alpha("white", 0.6), color = NA),
    legend.key = element_rect(fill = alpha("white", 0))
  ) +
  theme_natmed()

save_plot_safe(p_C_umap, "Figure1_C_umap", width = 4.8, height = 4)
