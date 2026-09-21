# ============================================================================
# Figure2_G_05_degree_threshold_sensitivity.R
# Project: Decoding Immune Landscapes and Therapeutic Targets in Intracerebral Hemorrhage
# Author:  Dimitrios Kyriakis
# ============================================================================

# ==============================================================================
# Figure 2, Panel G (internal check, not for manuscript) — is BHLHE41's low
# network degree (9, at the top-1% log-importance threshold) an artifact of
# that specific global threshold, or does it hold up across thresholds?
#
# The top-1% cutoff in Figure2_G_02_network_degree.R is applied GLOBALLY
# across all 35,515 motif-validated edges pooled together, not per-TF. That
# means it favors TFs whose individual edges happen to have high GRNBoost2
# importance scores -- which may correlate with the TF's own expression
# strength or motif prevalence -- over TFs with a genuinely broad but more
# moderately-weighted regulon. BHLHE41's FULL motif-validated regulon has
# 208 targets (not small), so this is a real, checkable question: does
# BHLHE41's rank improve substantially at looser thresholds, or is 9/208
# just where it genuinely sits?
#
# NOT a proposal to change the manuscript's threshold -- the top-1% choice
# already matches the S.Figure 4 legend and is defensible on its own terms.
# This is purely diagnostic: report what actually happens across a range of
# thresholds (including no threshold at all = full regulon size), and let
# the numbers say whatever they say. If BHLHE41 only looks good at a
# threshold chosen after seeing this result, that is not evidence for
# anything -- same logic as the PPARG significance-threshold discussion
# earlier in this project, just pointed at a TF we like this time.
#
# Input:  DERIVED_DIR/2026_08_ExtGenes_Regulons.csv        (pipelines/extract_loom_regulons.py)
#         DERIVED_DIR/2026_08_ExtGenes_adj.csv               (pipelines/scenic_GRN_ExtendedGenes.lsf.sh, GRNBoost2 step)
# Output: TABLE_DIR/Figure2_G_degree_threshold_sensitivity.csv   (every TF's out-degree AND rank at each threshold)
#         FIG_DIR/Figure2_G_degree_threshold_sensitivity.(pdf|png) (BHLHE41 and MAF's rank trajectory across thresholds)
# ==============================================================================

source("../../00_setup.R")
library(igraph)

THRESHOLDS <- c(1.00, 0.25, 0.10, 0.05, 0.01)  # fraction of edges KEPT (1.00 = full regulon, no cut; 0.01 = current top-1%)
HIGHLIGHT_TFS <- c("BHLHE41", "MAF", "JUNB", "FOS", "SPI1", "PPARG")

# ---- Validated (TF, target) edges + importance, same construction as
# Figure2_G_02_network_degree.R ----
regulons_df <- read.csv(file.path(DERIVED_DIR, "2026_08_ExtGenes_Regulons.csv"), row.names = 1, check.names = FALSE)
edge_idx <- which(as.matrix(regulons_df) == 1, arr.ind = TRUE)
validated_edges <- data.frame(
  TF     = gsub("[(+)]", "", colnames(regulons_df)[edge_idx[, "col"]]),
  target = rownames(regulons_df)[edge_idx[, "row"]],
  stringsAsFactors = FALSE
)

adj <- read.csv(file.path(DERIVED_DIR, "2026_08_ExtGenes_adj.csv"), check.names = FALSE)
edges <- merge(validated_edges, adj, by = c("TF", "target"))
edges$log_importance <- log1p(edges$importance)

# ---- Out-degree at each threshold (fraction of highest-importance edges kept) ----
results <- list()
for (frac in THRESHOLDS) {
  cutoff <- quantile(edges$log_importance, probs = 1 - frac)
  edges_kept <- edges[edges$log_importance >= cutoff, ]
  out_deg <- table(edges_kept$TF)
  df <- data.frame(TF = names(out_deg), out_degree = as.integer(out_deg), stringsAsFactors = FALSE)
  df$rank <- rank(-df$out_degree, ties.method = "min")
  df$threshold_pct <- frac * 100
  df$n_edges_total <- nrow(edges_kept)
  results[[as.character(frac)]] <- df
}
sensitivity <- bind_rows(results)

save_table_safe(sensitivity, "Figure2_G_degree_threshold_sensitivity.csv")

cat("---- Rank and out-degree of highlighted TFs across thresholds ----\n")
for (tf in HIGHLIGHT_TFS) {
  tf_rows <- sensitivity %>% filter(TF == tf) %>% arrange(desc(threshold_pct))
  if (nrow(tf_rows) == 0) {
    cat(tf, ": not present in the motif-validated edge set at all.\n")
    next
  }
  cat(sprintf("%-10s ", tf))
  for (i in seq_len(nrow(tf_rows))) {
    cat(sprintf("| top%2g%%: degree=%-4d rank=%-4d ", tf_rows$threshold_pct[i], tf_rows$out_degree[i], tf_rows$rank[i]))
  }
  cat("\n")
}

# ---- Plot: rank trajectory across thresholds for the highlighted TFs ----
plot_data <- sensitivity %>% filter(TF %in% HIGHLIGHT_TFS)
plot_data$threshold_pct <- factor(plot_data$threshold_pct, levels = sort(unique(plot_data$threshold_pct), decreasing = TRUE))

p_sensitivity <- ggplot(plot_data, aes(x = threshold_pct, y = rank, color = TF, group = TF)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  scale_y_reverse() +
  theme_bw() +
  labs(x = "% of highest-importance edges kept (100% = full regulon, no cut)",
       y = "Out-degree rank (1 = highest)",
       title = "Does BHLHE41's rank depend on the top-1% threshold choice?")

save_plot_safe(p_sensitivity, "Figure2_G_degree_threshold_sensitivity", width = 7, height = 5)

cat("\nSensitivity check complete. If BHLHE41's rank line stays flat/high across\n",
    "thresholds, its low top-1% rank is real, not a threshold artifact. If it\n",
    "climbs sharply as the threshold loosens, its full regulon IS broad, and the\n",
    "top-1% choice specifically is what's suppressing it.\n")
