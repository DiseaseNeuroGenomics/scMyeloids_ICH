# ============================================================================
# Figure2_G_FINAL_dotplot.R
# Project: Decoding Immune Landscapes and Therapeutic Targets in Intracerebral Hemorrhage
# Author:  Dimitrios Kyriakis
# ============================================================================

# ==============================================================================
# Figure 2, Panel G — FINAL recommended design. Combines everything this
# GRN validation effort established, without dropping anything that turned
# out to matter:
#
#   - Selection: differentially active TFs (p_val_adj < 0.05 on regulon
#     AUC), top 20 per direction by |log2FC| -- same population as v1/v4,
#     NOT additionally gated by Fisher significance. That extra gate is
#     what silently dropped PPARG earlier in this project even though it's
#     real (p_val_adj ~ 0) -- excluding it entirely would have hidden the
#     actual, reportable finding (differentially active, NOT enriched)
#     rather than shown it.
#   - x = TF, ranked by log2FC (ascending)
#   - y = log2FC
#   - color = log2FC (blue/lightgrey/red) -- matches the manuscript's
#     EXISTING Figure 2G legend convention, unchanged.
#   - size = Fisher OddsRatio -- also matches the existing legend ("size...
#     represents the odds ratio enrichment of the TF regulons in the
#     DEGs"), restored here and shown for every plotted TF, not just ones
#     that clear a significance gate -- so a low OddsRatio (PPARG, ~0.96)
#     reads as a small dot instead of being invisible.
#   - external validation (DoRothEA and/or CollecTri) shown as a thin tile
#     annotation strip above the dotplot (black = corroborated, grey =
#     novel), not as point shape -- standard genomics-figure convention
#     (same idea as a ComplexHeatmap top-annotation bar), aligned to the
#     same TF ordering as the dotplot below it via patchwork.
#
# The strict double-gated version (Figure2_G_TF_regulon_dotplot_ExtGenes.R)
# is untouched and still available for anyone who wants the more
# conservative population. This script is the recommended one for the
# manuscript's actual Figure 2G slot -- caption only needs one added
# sentence about the DoRothEA/CollecTri shape, everything else in the
# existing legend text still applies as written.
#
# Threshold note: uses p_val_adj < 0.05 (this project's locked FDR
# convention throughout), not the manuscript's original "adjusted p-value <
# 0.01" -- a deliberate choice, not an oversight; see the earlier
# discussion on fixing the raw-p / FDR handling across this whole panel.
#
# Input:  DERIVED_DIR/Figure2_G_DEG_TFs_ExtGenes.csv                        (Figure2_G_00_build_DEG_TFs_ExtGenes.R)
#         TABLE_DIR/Figure2_G_Fisher_TF_enrichment_ExtGenes_full.csv        (Figure2_G_TF_regulon_dotplot_ExtGenes.R)
#         TABLE_DIR/Figure2_G_dorothea_collectri_validation_TF_summary.csv  (Figure2_G_04_dorothea_collectri_validation.R)
# Output: TABLE_DIR/Figure2_G_FINAL_dotplot_data.csv
#         FIG_DIR/Figure2_G_FINAL_dotplot.(pdf|png)
# ==============================================================================

source("../../00_setup.R")
library(patchwork)

TOP_N_PER_SIDE <- 30

DEG_TFs <- read.csv(file.path(DERIVED_DIR, "Figure2_G_DEG_TFs_ExtGenes.csv"))
fisher_full <- read.csv(file.path(TABLE_DIR, "Figure2_G_Fisher_TF_enrichment_ExtGenes_full.csv"))
validation_summary <- read.csv(file.path(TABLE_DIR, "Figure2_G_dorothea_collectri_validation_TF_summary.csv"))

sig_TFs <- DEG_TFs %>% filter(p_val_adj < 0.05)

# ---- Merge in Fisher OddsRatio (unfiltered -- every tested TF, not just
# double-gated significant ones). TFs whose regulon never survived motif
# validation at all (absent from the Fisher table entirely) get OddsRatio =
# 1 (the null/no-enrichment default), so they still plot as a small dot
# rather than being silently dropped. ----
combined <- merge(sig_TFs, fisher_full[, c("TF_symbol", "OddsRatio", "FDR", "Significance")],
                   by.x = "gene", by.y = "TF_symbol", all.x = TRUE)
combined$OddsRatio[is.na(combined$OddsRatio)] <- 1
combined$Significance[is.na(combined$Significance)] <- "NS"

# ---- Merge in external validation status, graded by evidence strength, not
# collapsed to a flat yes/no. DoRothEA's A-E tiers reflect how many
# independent evidence lines back a given edge (A/B = curated literature +
# ChIP-seq + motif all agreeing; C = one strong source; D/E = motif-scan or
# expression-inference only, the weakest tier). A TF that only clears via a
# single weak CollecTri hit is not the same claim as one backed by DoRothEA
# tier A -- collapsing both to "Corroborated" would validate mere existence
# of prior evidence, not its quality. ----
combined <- merge(combined, validation_summary[, c("TF", "pct_in_dorothea", "pct_in_collectri", "best_dorothea_tier")],
                   by.x = "gene", by.y = "TF", all.x = TRUE)
combined$pct_in_dorothea[is.na(combined$pct_in_dorothea)] <- 0
combined$pct_in_collectri[is.na(combined$pct_in_collectri)] <- 0
combined$Evidence <- with(combined, ifelse(
  best_dorothea_tier %in% c("A", "B"), "DoRothEA A/B",
  ifelse(best_dorothea_tier == "C", "DoRothEA C",
  ifelse(best_dorothea_tier %in% c("D", "E"), "DoRothEA D/E",
  ifelse(pct_in_collectri > 0, "CollecTri only", "Novel")))))
combined$Evidence <- factor(combined$Evidence,
  levels = c("DoRothEA A/B", "DoRothEA C", "DoRothEA D/E", "CollecTri only", "Novel"))

save_table_safe(combined, "Figure2_G_FINAL_dotplot_data.csv")

# ---- Top N per direction, ranked by |standardized_diff_MTC456_minus_MTC123|
# (z-scored effect size), NOT raw |log2FC|. AUCell scores aren't calibrated
# across different regulons, so raw log2FC isn't directly comparable between
# TFs -- valid for describing one TF's own change, not for ranking TFs
# against each other. The standardized metric (from
# Figure2_G_00_build_DEG_TFs_ExtGenes.R) puts every TF on a common scale for
# exactly this selection step. log2FC is still what's plotted on the y-axis
# (more directly interpretable), just not what decides the top-N cut. ----
plot_data <- combined %>%
  mutate(direction = ifelse(log2FC > 0, "Up in MTC_456", "Up in MTC_123")) %>%
  group_by(direction) %>%
  slice_max(order_by = abs(standardized_diff_MTC456_minus_MTC123), n = TOP_N_PER_SIDE) %>%
  ungroup() %>%
  arrange(log2FC)
plot_data$Set <- factor(plot_data$gene, levels = unique(plot_data$gene))

# ---- Annotation strip: DoRothEA/CollecTri evidence strength as a thin tile
# bar above the dotplot (dark = strongest independent evidence, light =
# none), same x-order as plot_data$Set so it aligns exactly. ----
p_annot <- ggplot(plot_data, aes(x = Set, y = 1, fill = Evidence)) +
  geom_tile(color = "white", linewidth = 0.3) +
  scale_fill_manual(values = c(
    `DoRothEA A/B`   = "black",
    `DoRothEA C`     = "grey35",
    `DoRothEA D/E`   = "grey60",
    `CollecTri only` = "grey80",
    `Novel`          = "white"
  ), drop = FALSE) +
  theme_void() +
  theme(
    legend.position = "right",
    text = element_text(family = "Helvetica", size = 6)
  ) +
  labs(fill = "DoRothEA/CollecTri")

p_dots <- ggplot(plot_data, aes(x = Set, y = log2FC, size = OddsRatio, color = log2FC)) +
  geom_point() +
  theme_bw() +
  RotatedAxis() +
  scale_color_gradientn(
    colors = c("blue", "lightgrey", "red"),
    values = scales::rescale(c(min(plot_data$log2FC), 0, max(plot_data$log2FC)), midpoint = 0)
  ) +
  scale_size_continuous(range = c(0.3, 2.2)) +
  theme(
    text = element_text(family = "Helvetica", size = 6),
    legend.position = "right",
    legend.box = "horizontal",
    axis.text.x = element_text(angle = 90, hjust = 1)
  ) +
  guides(color = guide_colorbar(title.position = "top"), size = guide_legend(title.position = "top")) +
  labs(x = NULL, y = "log2FC (regulon AUC, MTC_456 vs MTC_123)", size = "Fisher OddsRatio")

# `&` (not `+`) applies the theme override to every sub-plot in the
# composition -- the standard patchwork way to control a collected guide's
# position, since guides="collect" alone doesn't expose a position argument.
p_final <- (p_annot / p_dots) +
  plot_layout(heights = c(1, 12), guides = "collect") &
  theme(legend.position = "bottom", legend.box = "horizontal",
        legend.key.size = unit(3, "mm"))

width_in  <- 170 * 0.03937
height_in <- 62  * 0.03937
save_plot_safe(p_final, "Figure2_G_FINAL_dotplot", width = width_in, height = height_in)

cat("Final Figure 2G:", nrow(plot_data), "TFs plotted (top", TOP_N_PER_SIDE, "per direction by |log2FC|).\n")
cat(" - Fisher-significant (FDR<0.05):", sum(plot_data$Significance == "Significant"), "\n")
cat(" - Evidence breakdown:\n")
print(table(plot_data$Evidence))
