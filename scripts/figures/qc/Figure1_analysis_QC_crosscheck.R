# ============================================================================
# Figure1_analysis_QC_crosscheck.R
# Project: Decoding Immune Landscapes and Therapeutic Targets in Intracerebral Hemorrhage
# Author:  Dimitrios Kyriakis
# ============================================================================

# ==============================================================================
# Figure 1 — QC re-derivation + cross-figure analysis (not a figure panel)
#
# Two things live here, both re-loading saved tables rather than recomputing
# from raw data, matching what the original notebook did at the end:
#   1. Sanity-check DEG counts / GSEA term counts against the numbers quoted
#      in the manuscript text.
#   2. Rebuild the population-level Myeloid GSEA and pull the leading-edge
#      genes of the suppressed lipid/lipoprotein biosynthesis terms — these
#      are consumed by Figure2_analysis_QC_crosscheck.R for the
#      MTC_456-overlap hypergeometric test.
#
# Input:  TABLE_DIR/Figure1_master_DEGs_all_assays.csv
#         TABLE_DIR/Figure1_F_myeloid_GSEA_full.csv
# Output: DERIVED_DIR/Figure1_leading_edge_list.rds  (consumed by Figure 2)
# ==============================================================================

source("../setup/00_setup.R")

DEGs_file <- read.csv(file.path(TABLE_DIR, "Figure1_master_DEGs_all_assays.csv"))
gsea_df   <- read.csv(file.path(TABLE_DIR, "Figure1_F_myeloid_GSEA_full.csv"))

# ---- 1. Sanity checks against manuscript-quoted numbers ----
cat("Myeloid DEGs (FDR<0.05):",
    DEGs_file %>% filter(adj.P.Val < 0.05, assay == "Myeloid") %>% nrow(), "\n")
cat("Unique significant genes across all assays (FDR<0.05):",
    length(unique(DEGs_file$ID[DEGs_file$adj.P.Val < 0.05])), "\n")

sig_counts <- DEGs_file %>%
  filter(adj.P.Val < 0.05) %>%
  count(assay, name = "n_sig_DEGs") %>%
  arrange(desc(n_sig_DEGs))
print(sig_counts)

sig_direction <- DEGs_file %>%
  filter(adj.P.Val < 0.05) %>%
  mutate(direction = ifelse(logFC > 0, "Unfavorable", "Favorable")) %>%
  count(assay, direction, name = "n") %>%
  tidyr::pivot_wider(names_from = direction, values_from = n, values_fill = 0)
print(sig_direction)

sig_terms <- gsea_df %>% filter(p.adjust < 0.05) %>% arrange(p.adjust)
cat("GO:BP terms enriched at FDR<0.05:", nrow(sig_terms), "\n")
save_table_safe(sig_terms, "Figure1_F_myeloid_GSEA_significant.csv")

# ---- 2. Rebuild population-level Myeloid GSEA, pull leading-edge genes ----
gmt_bp <- read.gmt(file.path(GMT_DIR, "c5.go.bp.v7.5.1.symbols.gmt"))

pop_data <- DEGs_file %>%
  filter(assay == "Myeloid", !is.na(t)) %>%
  arrange(desc(t))
gene_list <- setNames(pop_data$t, pop_data$ID)

set.seed(12345)
myeloid_res <- GSEA(
  geneList = gene_list, TERM2GENE = gmt_bp,
  minGSSize = 10, maxGSSize = 500, pvalueCutoff = 0.05, pAdjustMethod = "BH",
  eps = 0, seed = TRUE, nPermSimple = 10000
)
res_df <- as.data.frame(myeloid_res)

biosynthesis_terms <- c(
  "GOBP_LIPOPROTEIN_BIOSYNTHETIC_PROCESS",
  "GOBP_MEMBRANE_LIPID_BIOSYNTHETIC_PROCESS",
  "GOBP_GLYCOLIPID_BIOSYNTHETIC_PROCESS",
  "GOBP_MEMBRANE_LIPID_METABOLIC_PROCESS",
  "GOBP_LIPOSACCHARIDE_METABOLIC_PROCESS"
)

leading_edge_list <- res_df %>%
  filter(ID %in% biosynthesis_terms) %>%
  select(ID, core_enrichment) %>%
  { setNames(strsplit(.$core_enrichment, "/"), .$ID) }

print(leading_edge_list)

saveRDS(leading_edge_list, file.path(DERIVED_DIR, "Figure1_leading_edge_list.rds"))
cat("Saved: Figure1_leading_edge_list.rds (consumed by Figure2_analysis_QC_crosscheck.R)\n")

# ---- Reproducibility record ----
writeLines(capture.output(sessionInfo()), file.path(TABLE_DIR, "Figure1_sessionInfo.txt"))
cat("Saved: Figure1_sessionInfo.txt\n")
