# ============================================================================
# Figure2_G_00_build_DEG_TFs_ExtGenes.R
# Project: Decoding Immune Landscapes and Therapeutic Targets in Intracerebral Hemorrhage
# Author:  Dimitrios Kyriakis
# ============================================================================

# ==============================================================================
# Figure 2, Panel G (input build) — TF differential-activity table, ExtendedGenes
# panel: log2FC of each transcription factor's SCENIC regulon AUC score between
# MTC_456 and MTC_123.
#
# NEW script — does not touch Figure2_G_00_build_DEG_TFs.R (marked for
# deletion; reads the old pcaENTREZids-panel loom, 2025_MetaCL_auc_v10.loom,
# no longer the active gene panel). This version reads the ExtendedGenes
# panel's RegulonsAUC, extracted from 2026_08_ExtGenes_auc_v10.loom via
# pipelines/extract_loom_regulons.py (R's SCopeLoomR can't read that loom's
# compound HDF5 attributes directly — see that script's docstring).
#
# return.thresh = 1: FindAllMarkers's default (0.01) silently drops any TF
# with raw p_val >= 0.01 before this table is even written — a raw-p filter,
# not an FDR one. Kept off here so the FULL regulon table (every TF, with its
# p_val_adj) is what gets saved; significance filtering happens downstream,
# explicitly, on p_val_adj (see Figure2_G_TF_regulon_dotplot_ExtGenes.R).
# Matches the convention already used in Figure2_master_DEGs.R.
#
# Standardized effect size (added): AUCell scores are not calibrated across
# different regulons -- a regulon with a naturally wide AUC range can show a
# large raw log2FC without a proportionally large shift, and a
# tightly-distributed regulon can show a small raw log2FC despite a large
# shift relative to its own typical variation. That's invisible when
# looking at one TF alone, but this analysis repeatedly RANKS TFs against
# each other (top N by log2FC, "the largest shift of any TF") -- a
# comparison that needs scores on a common scale to be valid.
# `standardized_diff_MTC456_minus_MTC123` is that common-scale metric: each
# regulon's AUC is z-scored across ALL metacells (not just the two-cluster
# subset, so the scaling reference isn't circularly tied to the specific
# comparison being tested), then the mean z-score difference between groups
# is taken directly. This is NOT done by feeding z-scored values through
# Seurat's FindAllMarkers -- its avg_log2FC formula assumes log-normalized,
# non-negative count-like data (it calls expm1() to undo a log transform),
# which silently produces a meaningless number for arbitrary mean-zero
# input like z-scores. Computed here with plain rowMeans() instead.
# log2FC (raw AUC units) is kept unchanged and remains the right metric for
# describing a single TF's own magnitude of change; use
# standardized_diff_MTC456_minus_MTC123 for any claim comparing TFs to each
# other. Significance (p_val_adj) is unaffected by any of this -- Wilcoxon's
# rank-sum test is invariant to a monotonic per-feature transform.
#
# Input:  DERIVED_DIR/Figure2_Metacells.rds              (Figure2_A_umap_metacells.R)
#         DERIVED_DIR/2026_08_ExtGenes_RegulonsAUC.csv    (pipelines/extract_loom_regulons.py)
# Output: DERIVED_DIR/Figure2_G_DEG_TFs_ExtGenes.csv
# ==============================================================================

source("../../00_setup.R")

output_path <- file.path(DERIVED_DIR, "Figure2_G_DEG_TFs_ExtGenes.csv")

Metacells <- readRDS(file.path(DERIVED_DIR, "Figure2_Metacells.rds"))

# RegulonsAUC CSV is (cells x regulons) as written by pandas; Seurat assay
# data needs (features x cells), so transpose.
auc_df <- read.csv(file.path(DERIVED_DIR, "2026_08_ExtGenes_RegulonsAUC.csv"), row.names = 1, check.names = FALSE)
AUCmat <- t(as.matrix(auc_df))
rownames(AUCmat) <- gsub("[(+)]", "", rownames(AUCmat))

Metacells[["AUC"]] <- CreateAssayObject(data = AUCmat)

# Z-scored per regulon (row), across ALL metacells -- computed before
# subsetting to the two-cluster comparison, see rationale above.
AUC_z <- t(scale(t(AUCmat)))

subset_two_cl <- subset(Metacells, TwoClusters %in% c("MTC_123", "MTC_456"))
DefaultAssay(subset_two_cl) <- "AUC"
Idents(subset_two_cl) <- "TwoClusters"

DEG_TFs <- FindAllMarkers(subset_two_cl, only.pos = TRUE, logfc.threshold = 0, min.pct = 0.1,
                           return.thresh = 1, assay = "AUC")
DEG_TFs$TF_symbol <- DEG_TFs$gene
DEG_TFs$log2FC <- ifelse(DEG_TFs$cluster == "MTC_456", DEG_TFs$avg_log2FC, -1 * DEG_TFs$avg_log2FC)

# ---- Standardized mean difference, computed directly (not via Seurat) ----
cell_ids_456 <- colnames(subset_two_cl)[subset_two_cl$TwoClusters == "MTC_456"]
cell_ids_123 <- colnames(subset_two_cl)[subset_two_cl$TwoClusters == "MTC_123"]
mean_z_456 <- rowMeans(AUC_z[, cell_ids_456, drop = FALSE])
mean_z_123 <- rowMeans(AUC_z[, cell_ids_123, drop = FALSE])
standardized_effect <- data.frame(
  TF_symbol = names(mean_z_456),
  standardized_diff_MTC456_minus_MTC123 = mean_z_456 - mean_z_123[names(mean_z_456)],
  stringsAsFactors = FALSE
)

DEG_TFs <- merge(DEG_TFs, standardized_effect, by = "TF_symbol", all.x = TRUE)

write.csv(DEG_TFs, output_path, row.names = FALSE)
cat("Saved:", output_path, "(", nrow(DEG_TFs), "rows, all regulons, unfiltered )\n")
cat("Columns: log2FC (raw AUC units, for describing one TF's own change) and\n",
    "standardized_diff_MTC456_minus_MTC123 (z-scored, valid for ranking TFs against each other).\n")
