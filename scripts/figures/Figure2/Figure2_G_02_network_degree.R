# ============================================================================
# Figure2_G_02_network_degree.R
# Project: Decoding Immune Landscapes and Therapeutic Targets in Intracerebral Hemorrhage
# Author:  Dimitrios Kyriakis
# ============================================================================

# ==============================================================================
# Figure 2, Panel G (supplementary) — GRN network analysis, ExtendedGenes
# panel: rank TFs by weighted network centrality (motif-validated regulatory
# connections, weighted by GRNBoost2 importance) and by differential regulon
# activity (log2FC), and extract direct-neighbor subgraphs for a fixed set
# of TFs of interest: PPARG, SPI1, BHLHE41, TCF4.
#
# METHODOLOGY CHANGE from earlier versions of this script: previously used
# out-degree/betweenness on a network thresholded to the top 1% of
# log-transformed GRNBoost2 importance (matching the manuscript's original
# S.Figure 4 description). A sensitivity check
# (Figure2_G_05_degree_threshold_sensitivity.R) showed this specific global
# threshold systematically favors TFs with a few extreme-importance edges
# (JUNB: rank 1 only at the top-1% cut, mid-pack at every looser setting)
# over TFs with many solid, moderately-important edges (BHLHE41: rank 55
# with no threshold, rank 3 at top-10%, rank 16 at top-1% -- a real,
# threshold-dependent artifact). SPI1 -- the single largest log2FC hit in
# this analysis -- had ZERO edges surviving the top-1% cut and was entirely
# absent from the network centrality analysis under the old methodology.
#
# Corrected methodology:
#   - GRNBoost2 importance is a random-forest feature-importance score, not
#     a calibrated statistic (Moerman et al. 2019, Bioinformatics) -- no
#     principled statistical basis for an additional importance threshold
#     on top of motif validation.
#   - The statistically-principled edge filter in SCENIC is the
#     motif-validation (ctx) step itself, NES-based (Aibar et al. 2017,
#     Nature Methods) -- already applied upstream. Once an edge survives
#     that, it is the defensible network; no further arbitrary percentile
#     cut is applied here.
#   - Centrality is computed as WEIGHTED degree ("strength": sum of
#     GRNBoost2 importance across a TF's validated targets) and weighted
#     betweenness, the standard generalization of degree/betweenness to
#     weighted networks (Barrat, Barthelemy, Pastor-Satorras & Vespignani
#     2004, PNAS; Opsahl, Agneessens & Skvoretz 2010, Social Networks) --
#     avoids the threshold-artifact bias found above entirely, since nothing
#     is binarized.
#   - Current community best-practice guidance (sc-best-practices.org, GRN
#     chapter) explicitly recommends against arbitrary fixed-percentile
#     thresholds in GRN edge pruning.
#
# A top-N-by-importance cap is still used, but ONLY for the individual
# subgraph plots below (SUBGRAPH_TOP_N_EDGES), purely for display
# readability -- never for the quantitative ranking. That distinction was
# already correct in this script and is unchanged.
#
# Edge set = motif-validated (TF, target gene) pairs from
# 2026_08_ExtGenes_Regulons.csv, same file used throughout this panel -- NOT
# the raw, unvalidated GRNBoost2 adjacency matrix.
#
# Edge weight = GRNBoost2 importance, pulled back in from
# 2026_08_ExtGenes_adj.csv (TF, target, importance) for exactly the
# validated (TF, target) pairs.
#
# Input:  DERIVED_DIR/2026_08_ExtGenes_Regulons.csv        (pipelines/extract_loom_regulons.py)
#         DERIVED_DIR/2026_08_ExtGenes_adj.csv               (pipelines/scenic_GRN_ExtendedGenes.lsf.sh, GRNBoost2 step)
#         DERIVED_DIR/Figure2_G_DEG_TFs_ExtGenes.csv         (Figure2_G_00_build_DEG_TFs_ExtGenes.R)
# Output: TABLE_DIR/Figure2_G_network_degree_ranking.csv   (out_strength + weighted betweenness per TF)
#         FIG_DIR/Figure2_G_network_degree_barplot.(pdf|png)
#         FIG_DIR/Figure2_G_network_betweenness_barplot.(pdf|png)
#         FIG_DIR/Figure2_G_network_degree_barplot_updown_top10.(pdf|png)
#         FIG_DIR/Figure2_G_network_degree_barplot_updown_top20.(pdf|png)
#         FIG_DIR/Figure2_G_network_degree_barplot_updown_full.(pdf|png)
#         FIG_DIR/Figure2_G_network_degree_vs_betweenness.(pdf|png)
#         FIG_DIR/Figure2_G_network_degree_vs_log2FC.(pdf|png)
#         FIG_DIR/Figure2_G_subgraph_PPARG.(pdf|png)
#         FIG_DIR/Figure2_G_subgraph_SPI1.(pdf|png)
#         FIG_DIR/Figure2_G_subgraph_BHLHE41.(pdf|png)
#         FIG_DIR/Figure2_G_subgraph_TCF4.(pdf|png)
# ==============================================================================

source("../../setup/00_setup.R")
library(igraph)
library(ggraph)
library(ggrepel)

SUBGRAPH_TFS <- c("PPARG", "SPI1", "BHLHE41", "TCF4")
SUBGRAPH_TOP_N_EDGES <- 25

# ---- Validated (TF, target) edge set, from the consolidated regulons ----
regulons_df <- read.csv(file.path(DERIVED_DIR, "2026_08_ExtGenes_Regulons.csv"), row.names = 1, check.names = FALSE)
edge_idx <- which(as.matrix(regulons_df) == 1, arr.ind = TRUE)
validated_edges <- data.frame(
  TF     = gsub("[(+)]", "", colnames(regulons_df)[edge_idx[, "col"]]),
  target = rownames(regulons_df)[edge_idx[, "row"]],
  stringsAsFactors = FALSE
)
cat("Validated (motif-checked) TF-target edges:", nrow(validated_edges), "\n")

# ---- Pull GRNBoost2 importance back in for exactly these validated edges ----
adj <- read.csv(file.path(DERIVED_DIR, "2026_08_ExtGenes_adj.csv"), check.names = FALSE)
expected_cols <- c("TF", "target", "importance")
if (!all(expected_cols %in% colnames(adj))) {
  stop("2026_08_ExtGenes_adj.csv columns are: ", paste(colnames(adj), collapse = ", "),
       " -- expected: ", paste(expected_cols, collapse = ", "),
       ". Check the GRNBoost2/arboreto output schema before continuing.")
}

edges <- merge(validated_edges, adj, by = c("TF", "target"))
n_unmatched <- nrow(validated_edges) - nrow(edges)
if (n_unmatched > 0) {
  warning(n_unmatched, " validated edges had no matching importance score in adj.csv -- dropped.")
}

# ---- Build the FULL motif-validated graph -- no additional importance cut ----
g <- graph_from_data_frame(edges[, c("TF", "target", "importance")], directed = TRUE)
cat("Network:", vcount(g), "nodes,", ecount(g), "edges (full motif-validated set, no percentile cut).\n")

# ---- Weighted degree (strength): sum of GRNBoost2 importance across each
# TF's validated targets (Barrat et al. 2004) ----
tf_nodes <- unique(edges$TF)
out_str <- strength(g, vids = tf_nodes, mode = "out", weights = E(g)$importance)
degree_df <- data.frame(TF_symbol = names(out_str), out_strength = as.numeric(out_str), stringsAsFactors = FALSE) %>%
  arrange(desc(out_strength))

# ---- Weighted betweenness: edge weights as distances -- higher GRNBoost2
# importance means a shorter/easier path, so distance = 1/importance ----
bw <- betweenness(g, directed = TRUE, weights = 1 / E(g)$importance, normalized = FALSE)
degree_df$betweenness <- as.numeric(bw[degree_df$TF_symbol])
degree_df$betweenness[is.na(degree_df$betweenness)] <- 0

# ---- Merge with differential-activity log2FC ----
DEG_TFs <- read.csv(file.path(DERIVED_DIR, "Figure2_G_DEG_TFs_ExtGenes.csv"))
degree_df <- merge(degree_df, DEG_TFs[, c("gene", "log2FC", "p_val_adj")],
                    by.x = "TF_symbol", by.y = "gene", all.x = TRUE)
degree_df <- degree_df[order(-degree_df$out_strength), ]

save_table_safe(degree_df, "Figure2_G_network_degree_ranking.csv")

cat("---- Focal TFs: weighted degree (strength) and betweenness ----\n")
for (tf in SUBGRAPH_TFS) {
  tf_row <- degree_df[degree_df$TF_symbol == tf, ]
  if (nrow(tf_row) == 0) {
    cat(tf, ": not present in the motif-validated network at all.\n")
  } else {
    cat(sprintf("%-10s out_strength=%-12.2f betweenness=%-12.2f log2FC=%s\n",
                tf, tf_row$out_strength, tf_row$betweenness, round(tf_row$log2FC, 4)))
  }
}

# ---- Barplot: top TFs by out-strength ----
top_n_bar <- 25
degree_df_top <- head(degree_df, top_n_bar)
degree_df_top$TF_symbol <- factor(degree_df_top$TF_symbol, levels = rev(degree_df_top$TF_symbol))

p_degree_bar <- ggplot(degree_df_top, aes(x = TF_symbol, y = out_strength)) +
  geom_col(fill = "steelblue") +
  coord_flip() +
  theme_bw() +
  labs(x = "TF", y = "Out-strength (sum of GRNBoost2 importance, motif-validated targets)",
       title = paste("Top", top_n_bar, "TFs by weighted network degree (strength)"))

save_plot_safe(p_degree_bar, "Figure2_G_network_degree_barplot", width = 6, height = 8)

# ---- Barplot: top TFs by betweenness centrality ----
betweenness_df_top <- degree_df %>% arrange(desc(betweenness)) %>% head(top_n_bar)
betweenness_df_top$TF_symbol <- factor(betweenness_df_top$TF_symbol, levels = rev(betweenness_df_top$TF_symbol))

p_betweenness_bar <- ggplot(betweenness_df_top, aes(x = TF_symbol, y = betweenness)) +
  geom_col(fill = "darkorange") +
  coord_flip() +
  theme_bw() +
  labs(x = "TF", y = "Weighted betweenness centrality (motif-validated network)",
       title = paste("Top", top_n_bar, "TFs by betweenness centrality"))

save_plot_safe(p_betweenness_bar, "Figure2_G_network_betweenness_barplot", width = 6, height = 8)

# ---- Degree vs betweenness scatter ----
p_degree_vs_betweenness <- ggplot(degree_df, aes(x = out_strength, y = betweenness)) +
  geom_point(aes(color = TF_symbol %in% SUBGRAPH_TFS), size = 1.5, alpha = 0.7) +
  scale_color_manual(values = c(`FALSE` = "grey70", `TRUE` = "firebrick"), guide = "none") +
  geom_text_repel(data = subset(degree_df, TF_symbol %in% SUBGRAPH_TFS), aes(label = TF_symbol), size = 3) +
  theme_bw() +
  labs(x = "Out-strength (weighted degree)", y = "Betweenness centrality (bridging role)")

save_plot_safe(p_degree_vs_betweenness, "Figure2_G_network_degree_vs_betweenness", width = 6, height = 5)

# ---- Barplot: top TFs by out-strength, split by direction of differential
# activity (MTC_456 vs MTC_123) -- three flavors ----
degree_updown_base <- degree_df %>%
  filter(!is.na(log2FC)) %>%
  mutate(direction = ifelse(log2FC > 0, "Up in MTC_456", "Up in MTC_123"))

top_n_per_side_degree <- function(df, n) {
  df %>%
    group_by(direction) %>%
    slice_max(order_by = out_strength, n = n) %>%
    ungroup()
}

make_degree_barplot_updown <- function(df, out_name, title) {
  df <- df %>%
    group_by(direction) %>%
    mutate(TF_symbol = factor(TF_symbol, levels = TF_symbol[order(out_strength)])) %>%
    ungroup()

  p <- ggplot(df, aes(x = TF_symbol, y = out_strength, fill = direction)) +
    geom_col() +
    coord_flip() +
    facet_wrap(~direction, scales = "free_y") +
    scale_fill_manual(values = c(`Up in MTC_456` = "firebrick", `Up in MTC_123` = "steelblue"), guide = "none") +
    theme_bw() +
    labs(x = "TF", y = "Out-strength (sum of GRNBoost2 importance)", title = title)

  height_in <- max(4, max(table(df$direction)) * 0.3)
  save_plot_safe(p, out_name, width = 9, height = height_in)
  cat(out_name, ":", nrow(df), "TFs plotted.\n")
}

make_degree_barplot_updown(top_n_per_side_degree(degree_updown_base, 10), "Figure2_G_network_degree_barplot_updown_top10",
                            "Top 10 TFs by out-strength per direction (MTC_456 vs MTC_123)")

make_degree_barplot_updown(top_n_per_side_degree(degree_updown_base, 20), "Figure2_G_network_degree_barplot_updown_top20",
                            "Top 20 TFs by out-strength per direction (MTC_456 vs MTC_123)")

make_degree_barplot_updown(degree_updown_base, "Figure2_G_network_degree_barplot_updown_full",
                            "All network TFs by out-strength per direction (MTC_456 vs MTC_123)")

# ---- Degree vs log2FC scatter, highlighting the same 4 TFs of interest ----
degree_df$highlight <- degree_df$TF_symbol %in% SUBGRAPH_TFS

p_degree_fc <- ggplot(degree_df, aes(x = out_strength, y = log2FC)) +
  geom_point(aes(color = highlight), size = 1.5, alpha = 0.7) +
  scale_color_manual(values = c(`FALSE` = "grey70", `TRUE` = "firebrick"), guide = "none") +
  geom_text_repel(data = subset(degree_df, highlight), aes(label = TF_symbol), size = 3) +
  theme_bw() +
  labs(x = "Out-strength (weighted network degree)", y = "log2FC (regulon AUC, MTC_456 vs MTC_123)")

save_plot_safe(p_degree_fc, "Figure2_G_network_degree_vs_log2FC", width = 6, height = 5)

# ---- Direct-neighbor subgraphs for PPARG, SPI1, BHLHE41, TCF4 ----
# Uses the full validated edge set (unchanged by this revision) -- capped
# per-TF to top N targets by importance purely for display readability,
# consistent with "cap for visualization only, never for the quantitative
# ranking," applied consistently across this whole script now.
for (tf in SUBGRAPH_TFS) {
  tf_edges <- edges %>% filter(TF == tf) %>% slice_max(order_by = importance, n = SUBGRAPH_TOP_N_EDGES)
  if (nrow(tf_edges) == 0) {
    warning(tf, " has no motif-validated regulon at all (absent from 2026_08_ExtGenes_Regulons.csv) -- skipping subgraph.")
    next
  }
  cat(tf, ": plotting", nrow(tf_edges), "of", sum(edges$TF == tf), "total validated targets.\n")

  sub_g <- graph_from_data_frame(tf_edges[, c("TF", "target", "importance")], directed = TRUE)

  p_sub <- ggraph(sub_g, layout = "star", center = which(V(sub_g)$name == tf)) +
    geom_edge_link(aes(width = importance), alpha = 0.6, color = "grey30", arrow = arrow(length = unit(2, "mm"))) +
    geom_node_point(aes(color = name == tf), size = 4) +
    geom_node_text(aes(label = name), repel = TRUE, size = 3, color = "black") +
    scale_color_manual(values = c(`FALSE` = "grey40", `TRUE` = "firebrick"), guide = "none") +
    theme_void() +
    theme(plot.background = element_rect(fill = "white", color = NA),
          panel.background = element_rect(fill = "white", color = NA)) +
    labs(title = paste(tf, "direct neighbors (top", nrow(tf_edges), "validated targets by importance)"))

  save_plot_safe(p_sub, paste0("Figure2_G_subgraph_", tf), width = 6, height = 6)
}

cat("Network analysis complete.\n")
