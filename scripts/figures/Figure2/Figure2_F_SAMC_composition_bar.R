# ============================================================================
# Figure2_F_SAMC_composition_bar.R
# Project: Decoding Immune Landscapes and Therapeutic Targets in Intracerebral Hemorrhage
# Author:  Dimitrios Kyriakis
# ============================================================================

# ==============================================================================
# Figure 2, Panel F — mapping human MTC metacells onto the mouse stroke
# reference (Beuker et al. 2022): stacked bar of predicted mouse cell-type
# composition per MTC cluster, testing whether mouse stroke-associated
# myeloid cells (SAMC) align preferentially with MTC_456.
#
# Input:  DERIVED_DIR/Figure2_F_SAMC_reference.rds (Figure2_F_00_build_SAMC_reference.R)
#         DERIVED_DIR/Figure2_Metacells.rds        (Figure2_A_umap_metacells.R)
# Output: TABLE_DIR/Figure2_F_predicted_celltype_composition.csv
#         FIG_DIR/Figure2_F_SAMC_composition_bar.(pdf|png)
#         FIG_DIR/Figure2_F_SAMC_reference_and_projection.(pdf|png)  (context panel)
# ==============================================================================

source("../../setup/00_setup.R")
library(dittoSeq)

colors_SAMC <- c(
  "#FFC312", "#C4E538", "#12CBC4", "#FDA7DF", "#ED4C67", "#F79F1F", "#A3CB38",
  "#1289A7", "#D980FA", "#B53471", "#EE5A24", "#009432", "#0652DD", "#9980FA",
  "#833471", "#EA2027", "#006266", "#1B1464", "#5758BB", "#6F1E51", "#40407A"
)
names(colors_SAMC) <- c(
  "Micro_1", "Micro_2", "Micro_3", "stress_Micro", "CAM_1", "CAM_2", "SAMC",
  "Macro_1", "Macro_2", "stress_Myeloid", "mDC1", "mDC2", "Granulo_1", "Granulo_2",
  "Mast", "prolif_cells", "Bc", "gdTc", "ILC2", "Tc", "NK"
)

ifnb      <- readRDS(file.path(DERIVED_DIR, "Figure2_F_SAMC_reference.rds"))
Metacells <- readRDS(file.path(DERIVED_DIR, "Figure2_Metacells.rds"))

anchors <- FindTransferAnchors(
  reference = ifnb, query = Metacells, normalization.method = "SCT",
  reference.reduction = "pca", dims = 1:40
)
Metacells <- MapQuery(
  anchorset = anchors, query = Metacells, reference = ifnb,
  refdata = list(celltype = "cluster"), reference.reduction = "harmony",
  reduction.model = "umap.harmony"
)

p_ref     <- DimPlot(ifnb, group.by = "cluster", reduction = "umap.harmony", pt.size = 0.0001,
                      cols = colors_SAMC, label = FALSE) + xlim(-15, 10) + ylim(-15, 10) + ggtitle("Reference")
p_project <- DimPlot(Metacells, reduction = "umap", group.by = "predicted.celltype", cols = colors_SAMC)
p_bar_ctx <- dittoBarPlot(Metacells, group.by = "Clusters", var = "predicted.celltype",
                           color.panel = colors_SAMC) + coord_flip()

p_context <- cowplot::plot_grid(p_ref, p_project, p_bar_ctx, align = "h", ncol = 3, rel_widths = c(1, 1, 1))
save_plot_safe(p_context, "Figure2_F_SAMC_reference_and_projection", width = 18, height = 5)

# ---- Panel F: predicted mouse cell-type composition per MTC cluster ----
p_F_bar <- dittoBarPlot(Metacells, "predicted.celltype", group.by = "Clusters", color.panel = colors_SAMC) +
  theme(text = element_text(size = 20),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

save_plot_safe(p_F_bar, "Figure2_F_SAMC_composition_bar", width = 6.6, height = 5)

composition_table <- Metacells@meta.data %>%
  count(Clusters, predicted.celltype) %>%
  group_by(Clusters) %>%
  mutate(fraction = n / sum(n)) %>%
  ungroup()
save_table_safe(composition_table, "Figure2_F_predicted_celltype_composition.csv")
