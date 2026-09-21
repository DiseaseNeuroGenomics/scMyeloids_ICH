# ============================================================================
# Figure2_A_umap_metacells.R
# Project: Decoding Immune Landscapes and Therapeutic Targets in Intracerebral Hemorrhage
# Author:  Dimitrios Kyriakis
# ============================================================================

# ==============================================================================
# Figure 2, Panel A — Myeloid metacell UMAP (scANVI transfer-annotation +
# hdWGCNA cell aggregation), colored by MTC subclass
#
# Input:  BASE_DIR/2024_03_28_Myeloid_Metacells_Subclass_ADAM.rds.ztsd
# Output: FIG_DIR/Figure2_A_umap_all_clusters.(pdf|png)
#         FIG_DIR/Figure2_A_umap_clusters_labeled.(pdf|png)
#         DERIVED_DIR/Figure2_Metacells.rds  (Metacells + TwoClusters column,
#           reused by every other Figure2_*.R script so the ~large Seurat
#           object and cluster labels aren't re-derived each time)
# ==============================================================================

source("../../00_setup.R")
source("../../00_helpers_stats.R")
source("../../00_theme_colors.R")

input_rds <- file.path(BASE_DIR, "2024_03_28_Myeloid_Metacells_Subclass_ADAM.rds.ztsd")

Metacells <- readCRDS(input_rds)

Metacells$Clusters <- factor(
  Metacells$Clusters,
  levels = c("MTC_1", "MTC_2", "MTC_3", "MTC_4", "MTC_5", "MTC_6", "MTC_Prolif")
)

Metacells$TwoClusters <- as.vector(Metacells$Clusters)
Metacells$TwoClusters[Metacells$Clusters %in% c("MTC_1", "MTC_2", "MTC_3")] <- "MTC_123"
Metacells$TwoClusters[Metacells$Clusters %in% c("MTC_4", "MTC_5", "MTC_6")] <- "MTC_456"

p_A_umap_full <- DimPlot(Metacells, group.by = "Clusters", cols = MTC_Colors) +
  theme_natmed()
p_A_umap_labeled <- DimPlot(Metacells, group.by = "Clusters", label = TRUE, raster = TRUE) +
  theme_natmed()

save_plot_safe(p_A_umap_full,    "Figure2_A_umap_all_clusters",    width = 6, height = 5)
save_plot_safe(p_A_umap_labeled, "Figure2_A_umap_clusters_labeled", width = 6, height = 5)

saveRDS(Metacells, file.path(DERIVED_DIR, "Figure2_Metacells.rds"))
cat("Saved: Figure2_Metacells.rds\n")
