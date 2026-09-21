# ============================================================================
# Figure2_G_04_dorothea_collectri_validation.R
# Project: Decoding Immune Landscapes and Therapeutic Targets in Intracerebral Hemorrhage
# Author:  Dimitrios Kyriakis
# ============================================================================

# ==============================================================================
# Figure 2, Panel G (validation) — cross-reference our GRNBoost2 + pySCENIC
# ctx motif-validated TF-target edges (ExtendedGenes panel) against two
# independent, literature/database-curated TF-target resources:
#
#   DoRothEA  — confidence-tiered A-E. A/B = curated literature databases +
#               ChIP-seq binding evidence + motif support agreeing; C =
#               single strong evidence line; D/E = motif-scan or
#               expression-inference only (weakest tier).
#   CollecTri — a more recent, broader aggregation of ~12 TF-target
#               resources (incl. DoRothEA, TRRUST, ExTRI). No A-E tiering,
#               but usually carries a per-edge list of supporting databases.
#
# Purpose: our own pipeline's validation (GRNBoost2 co-expression + pySCENIC
# motif enrichment) is entirely internal to this dataset -- no independent
# check against literature or ChIP-seq evidence. This script adds that
# external check. ANNOTATION only -- does not filter or change our own
# edges/TFs.
#
# Fetched as flat TSV directly from the OmniPath REST API (omnipathdb.org),
# NOT via the decoupleR/OmnipathR R packages -- OmnipathR failed to install
# (hard-requires xml2 >= 1.4.0, cluster has 1.3.6; forcing that upgrade
# risked rippling into other packages that depend on xml2, e.g. httr/rvest,
# which we deliberately avoided touching). This route needs zero new R
# packages, just outbound HTTPS from wherever this runs.
#
# Endpoints verified directly (2026-08-07) against the live API:
#   DoRothEA:  https://omnipathdb.org/interactions?datasets=dorothea&genesymbols=1&fields=dorothea_level&format=tsv
#              columns: source, target, source_genesymbol, target_genesymbol,
#              is_directed, is_stimulation, is_inhibition, consensus_direction,
#              consensus_stimulation, consensus_inhibition, dorothea_level
#   CollecTri: https://omnipathdb.org/interactions?datasets=collectri&genesymbols=1&fields=databases&format=tsv
#              same base columns; `databases` field requested but CONFIRMED
#              ABSENT from the actual response (checked the printed column
#              list from a real run) -- so per-edge resource-count isn't
#              available from this endpoint. Membership in CollecTri is
#              still correctly detected via an explicit flag set on every
#              row of the raw table, independent of that missing field.
#
# IF download.file() FAILS (no outbound internet from this node -- common on
# HPC compute nodes, less common on login nodes, unconfirmed for this
# cluster): download the two URLs above from your own machine (browser or
# curl) and scp them to DERIVED_DIR as dorothea_omnipath.tsv and
# collectri_omnipath.tsv, then rerun -- the script skips the download if the
# files already exist.
#
# Input:  DERIVED_DIR/2026_08_ExtGenes_Regulons.csv        (pipelines/extract_loom_regulons.py; our validated edges)
#         TABLE_DIR/Figure2_G_network_degree_ranking.csv     (Figure2_G_02_network_degree.R; degree/betweenness/log2FC per TF, for context)
#         DERIVED_DIR/dorothea_omnipath.tsv                  (downloaded here if absent)
#         DERIVED_DIR/collectri_omnipath.tsv                 (downloaded here if absent)
# Output: TABLE_DIR/Figure2_G_dorothea_collectri_validation_edges.csv       (every one of our edges, annotated)
#         TABLE_DIR/Figure2_G_dorothea_collectri_validation_TF_summary.csv (per-TF corroboration %, merged with degree/betweenness/log2FC)
#         FIG_DIR/Figure2_G_dorothea_confidence_distribution.(pdf|png)
# ==============================================================================

source("../../00_setup.R")

DOROTHEA_URL  <- "https://omnipathdb.org/interactions?datasets=dorothea&genesymbols=1&fields=dorothea_level&format=tsv"
COLLECTRI_URL <- "https://omnipathdb.org/interactions?datasets=collectri&genesymbols=1&fields=databases&format=tsv"

dorothea_path  <- file.path(DERIVED_DIR, "dorothea_omnipath.tsv")
collectri_path <- file.path(DERIVED_DIR, "collectri_omnipath.tsv")

if (!file.exists(dorothea_path)) {
  cat("Downloading DoRothEA interactions from OmniPath...\n")
  download.file(DOROTHEA_URL, dorothea_path, method = "curl")
}
if (!file.exists(collectri_path)) {
  cat("Downloading CollecTri interactions from OmniPath...\n")
  download.file(COLLECTRI_URL, collectri_path, method = "curl")
}

dorothea_net  <- read.delim(dorothea_path, stringsAsFactors = FALSE)
collectri_net <- read.delim(collectri_path, stringsAsFactors = FALSE)
cat("DoRothEA columns:", paste(colnames(dorothea_net), collapse = ", "), "\n")
cat("CollecTri columns:", paste(colnames(collectri_net), collapse = ", "), "\n")

# ---- Our validated (TF, target) edges, ExtendedGenes panel ----
regulons_df <- read.csv(file.path(DERIVED_DIR, "2026_08_ExtGenes_Regulons.csv"), row.names = 1, check.names = FALSE)
edge_idx <- which(as.matrix(regulons_df) == 1, arr.ind = TRUE)
our_edges <- data.frame(
  TF     = gsub("[(+)]", "", colnames(regulons_df)[edge_idx[, "col"]]),
  target = rownames(regulons_df)[edge_idx[, "row"]],
  stringsAsFactors = FALSE
)
cat("Our validated (GRNBoost2 + motif) TF-target edges:", nrow(our_edges), "\n")

# ---- DoRothEA lookup (gene-symbol columns, not the UniProt source/target IDs) ----
dorothea_lookup <- dorothea_net %>%
  transmute(TF = source_genesymbol, target = target_genesymbol, dorothea_confidence = dorothea_level) %>%
  distinct()

# ---- CollecTri lookup. `databases` field confirmed NOT present in this
# API response (checked the printed column list) -- existence in CollecTri
# is tracked via an explicit flag (in_collectri_table) set on every row of
# collectri_net itself, independent of that missing field, so presence is
# still correctly detected. collectri_n_resources stays NA (no per-edge
# resource-count data available from this endpoint) -- it's a bonus field,
# not what determines membership. ----
collectri_lookup <- collectri_net %>%
  transmute(TF = source_genesymbol, target = target_genesymbol,
            resources = if ("databases" %in% colnames(collectri_net)) databases else NA_character_) %>%
  distinct()
collectri_lookup$in_collectri_table <- TRUE
collectri_lookup$collectri_n_resources <- ifelse(
  is.na(collectri_lookup$resources), NA_integer_,
  lengths(strsplit(collectri_lookup$resources, "[;,]"))
)

# ---- Annotate every one of our edges ----
validated <- our_edges %>%
  left_join(dorothea_lookup, by = c("TF", "target")) %>%
  left_join(collectri_lookup %>% select(TF, target, in_collectri_table, collectri_n_resources), by = c("TF", "target"))

validated$in_dorothea  <- !is.na(validated$dorothea_confidence)
validated$in_collectri <- !is.na(validated$in_collectri_table)
validated$independently_corroborated <- validated$in_dorothea | validated$in_collectri

save_table_safe(validated, "Figure2_G_dorothea_collectri_validation_edges.csv")

cat("\nOf", nrow(validated), "of our validated edges:\n")
cat(" - in DoRothEA (any tier):", sum(validated$in_dorothea),
    sprintf("(%.1f%%)\n", 100 * mean(validated$in_dorothea)))
cat(" - in DoRothEA tier A/B:", sum(validated$dorothea_confidence %in% c("A", "B"), na.rm = TRUE), "\n")
cat(" - in CollecTri:", sum(validated$in_collectri),
    sprintf("(%.1f%%)\n", 100 * mean(validated$in_collectri)))
cat(" - either resource:", sum(validated$independently_corroborated),
    sprintf("(%.1f%%)\n", 100 * mean(validated$independently_corroborated)))

# ---- Per-TF summary, merged with degree/betweenness/log2FC for context ----
tf_summary <- validated %>%
  group_by(TF) %>%
  summarise(
    n_targets = n(),
    pct_in_dorothea = mean(in_dorothea) * 100,
    pct_in_collectri = mean(in_collectri) * 100,
    pct_corroborated = mean(independently_corroborated) * 100,
    best_dorothea_tier = if (any(in_dorothea)) min(dorothea_confidence[in_dorothea]) else NA_character_,
    .groups = "drop"
  ) %>%
  arrange(desc(pct_corroborated))

degree_ranking <- read.csv(file.path(TABLE_DIR, "Figure2_G_network_degree_ranking.csv"))
tf_summary <- merge(tf_summary, degree_ranking[, c("TF_symbol", "out_strength", "betweenness", "log2FC")],
                     by.x = "TF", by.y = "TF_symbol", all.x = TRUE)

save_table_safe(tf_summary, "Figure2_G_dorothea_collectri_validation_TF_summary.csv")

cat("\n---- Focal TFs (PPARG, SPI1, BHLHE41, TCF4) ----\n")
print(tf_summary %>% filter(TF %in% c("PPARG", "SPI1", "BHLHE41", "TCF4")))

# ---- Barplot: DoRothEA confidence-tier distribution among our edges.
# "Not in DoRothEA" is its own bar (not dropped), so the plot is honest
# about coverage, not just tier quality among the subset that happens to
# match. ----
tier_counts <- validated %>%
  mutate(tier = ifelse(in_dorothea, dorothea_confidence, "Not in DoRothEA")) %>%
  count(tier) %>%
  mutate(tier = factor(tier, levels = c("A", "B", "C", "D", "E", "Not in DoRothEA")))

p_tiers <- ggplot(tier_counts, aes(x = tier, y = n)) +
  geom_col(fill = "steelblue") +
  theme_bw() +
  labs(x = "DoRothEA confidence tier", y = "Number of our validated edges",
       title = "How many of our GRNBoost2 + motif-validated edges are\nindependently corroborated by DoRothEA?")

save_plot_safe(p_tiers, "Figure2_G_dorothea_confidence_distribution", width = 6, height = 4)

cat("\nValidation complete.\n")
