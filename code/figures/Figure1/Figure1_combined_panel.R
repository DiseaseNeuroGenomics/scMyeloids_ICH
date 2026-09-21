# ============================================================================
# Figure1_combined_panel.R
# Project: Decoding Immune Landscapes and Therapeutic Targets in Intracerebral Hemorrhage
# Author:  Dimitrios Kyriakis
# ============================================================================

# ==============================================================================
# Figure 1 — combined C+D+E+F row (final assembled panel, minus the schematic
# (A) and CCA metadata heatmap (B), which have no source code in this project)
#
# Re-generates each panel plot from its own script's saved objects on disk
# rather than re-running the analyses, then arranges them side by side.
#
# Input:  Figure1_C_umap.R, Figure1_D_mRS_effect_forest.R,
#         Figure1_E_deg_summary_bar.R, Figure1_F_myeloid_enrichment.R
#         (run once each beforehand; this script re-sources them for the plot
#         objects, which is cheap for C, and re-derives D/E/F from their
#         already-cached inputs)
# Output: FIG_DIR/Figure1_combined_CDEF.(pdf|png)
# ==============================================================================

source("../../00_setup.R")
source("../../00_helpers_stats.R")
init_natmed_fonts()

# Re-source each panel script to rebuild its plot object in this session.
# (Each script is independent and safe to source standalone; sourcing them
# here just avoids duplicating panel-construction code.)
source("Figure1_C_umap.R")               # -> p_C_umap
source("Figure1_D_mRS_effect_forest.R")  # -> p_D_forest
source("Figure1_E_deg_summary_bar.R")    # -> p_E_bar
source("Figure1_F_myeloid_enrichment.R") # -> p_F_gsea

combined <- wrap_plots(list(p_C_umap, p_D_forest, p_E_bar, p_F_gsea),
                        nrow = 1, widths = c(1.2, 1, 1, 1))

widths <- c(p_C_umap = 4.8, p_D_forest = 4, p_E_bar = 4, p_F_gsea = 4)
save_plot_safe(combined, "Figure1_combined_CDEF", width = sum(widths), height = 4)
