# ============================================================================
# Figure2_G_TF_regulon_dotplot_ExtGenes.R
# Project: Decoding Immune Landscapes and Therapeutic Targets in Intracerebral Hemorrhage
# Author:  Dimitrios Kyriakis
# ============================================================================

# ==============================================================================
# Figure 2, Panel G — differentially active transcription factors (TFs)
# between MTC_456 and MTC_123, ExtendedGenes panel, sized by Fisher's-exact
# enrichment of each TF's SCENIC regulon among the combined ("Disregulated")
# MTC_456/MTC_123 marker gene set.
#
# NEW script — does not touch Figure2_G_TF_regulon_dotplot.R (old
# pcaENTREZids panel, marked for deletion) or Figure2_G_ALT_ExtendedGenesLoom.R
# (diagnostic combined script this was split out of, also marked for
# deletion). Two fixes relative to both of those:
#
# 1. Significance is FDR-based throughout, not raw p_val:
#      - DEGs_genes (feeding the "Disregulated" gene set): p_val_adj < 0.05,
#        not p_val < 0.01.
#      - TF-level significance gate (sig_TFs): explicit p_val_adj < 0.05 on
#        Figure2_G_DEG_TFs_ExtGenes.csv, since that table is now unfiltered
#        (return.thresh = 1 upstream) rather than implicitly pre-filtered by
#        FindAllMarkers's raw-p default.
#    (The Fisher-test-derived "Significance" column was already FDR/BH-based
#    in both prior versions — unchanged here.)
#
# 2. Background gene universe fixed. Figure2_G_ALT_ExtendedGenesLoom.R
#    computed Background AFTER switching subset_two_cl's DefaultAssay to
#    "AUC", so rownames(subset_two_cl) returned the ~200 regulon names
#    instead of genes — an invalid Fisher-test background (regulon vs
#    regulon, not gene vs gene). Background here is instead the actual
#    ExtendedGenes SCENIC input gene universe (rownames of the Regulons
#    incidence matrix, ~6,246 genes) — the correct "genes that had a chance
#    to be nominated as a regulon target" universe.
#
# Input:  DERIVED_DIR/Figure2_G_DEG_TFs_ExtGenes.csv        (Figure2_G_00_build_DEG_TFs_ExtGenes.R)
#         DERIVED_DIR/2026_08_ExtGenes_Regulons.csv          (pipelines/extract_loom_regulons.py)
#         TABLE_DIR/Figure2_master_DEGs_FindAllMarkers.csv   (Figure2_master_DEGs.R)
# Output: TABLE_DIR/Figure2_G_Fisher_TF_enrichment_ExtGenes_full.csv
#         TABLE_DIR/Figure2_G_Fisher_TF_enrichment_ExtGenes_significant.csv
#         FIG_DIR/Figure2_G_TF_regulon_dotplot_ExtGenes.(pdf|png)
# ==============================================================================

source("../../setup/00_setup.R")

# regulonsToGeneLists(): inlined from aertslab/SCENIC's R/aux_regulons.R —
# see Figure2_G_TF_regulon_dotplot.R for the full rationale. incidMat here:
# TFs as rows, genes as columns.
regulonsToGeneLists <- function(incidMat) {
  regulons <- list()
  for (i in rownames(incidMat))
    regulons[[i]] <- colnames(incidMat)[which(incidMat[i, ] == 1)]
  return(regulons)
}

# ---- Regulon gene lists + gene universe (Background) ----
# CSV is (genes x regulons); regulonsToGeneLists() expects (regulons x genes),
# so transpose for TF_list. Background is genes (rows), read BEFORE transpose.
regulons_df <- read.csv(file.path(DERIVED_DIR, "2026_08_ExtGenes_Regulons.csv"), row.names = 1, check.names = FALSE)
Background <- rownames(regulons_df)
cat("Fisher-test background gene universe:", length(Background), "genes\n")
regulons_incidMat <- t(as.matrix(regulons_df))
TF_list <- regulonsToGeneLists(regulons_incidMat)

# ---- DEG_TFs (unfiltered, from build script) + FDR-based TF significance gate ----
DEG_TFs <- read.csv(file.path(DERIVED_DIR, "Figure2_G_DEG_TFs_ExtGenes.csv"))
sig_TFs <- DEG_TFs$gene[DEG_TFs$p_val_adj < 0.05]

# ---- DEGs_genes: FDR-based (p_val_adj < 0.05), not raw p_val < 0.01 ----
DEGs_genes <- read.csv(file.path(TABLE_DIR, "Figure2_master_DEGs_FindAllMarkers.csv")) %>%
  filter(p_val_adj < 0.05)

up_genes   <- DEGs_genes$gene[DEGs_genes$cluster == "MTC_456"]
down_genes <- DEGs_genes$gene[DEGs_genes$cluster == "MTC_123"]
DEGs_list  <- list("Disregulated" = unique(c(up_genes, down_genes)))

# ------------------------------------------------------------------------------
# Fisher's exact test: is each TF's regulon enriched among the combined DEG set?
# ------------------------------------------------------------------------------
Calculate_Fisher_Enrichement2 <- function(List1, List2, Background) {
  A <- length(intersect(List1, List2))
  B <- length(List1[!List1 %in% List2])
  C <- length(List2[!List2 %in% List1])
  D <- length(Background[!Background %in% c(List1, List2)])
  observed_data <- matrix(c(A, B, C, D), nrow = 2)
  fisher_result <- fisher.test(observed_data)
  list(p_value = fisher_result$p.value, odd_ratio = as.numeric(fisher_result$estimate))
}

perform_fisher_tests <- function(DEGs_list, TF_list, background_genes, fdr_threshold = 0.05) {
  results <- data.frame(
    Reference = character(), TF_name = character(), TF_symbol = character(),
    GeneRatio = numeric(), OddsRatio = numeric(), PValue = numeric(),
    FDR = numeric(), Significance = character(), stringsAsFactors = FALSE
  )

  for (deg_name in names(DEGs_list)) {
    deg_genes <- DEGs_list[[deg_name]]
    for (tf_name in names(TF_list)) {
      tf_genes <- TF_list[[tf_name]]
      tf_symbol <- gsub("[(+)]", "", tf_name)
      fisher_result <- Calculate_Fisher_Enrichement2(deg_genes, tf_genes, background_genes)
      gene_ratio <- length(intersect(deg_genes, tf_genes)) / length(deg_genes)
      results <- rbind(results, data.frame(
        Reference = deg_name, TF_name = tf_name, TF_symbol = tf_symbol,
        GeneRatio = gene_ratio, OddsRatio = fisher_result$odd_ratio, PValue = fisher_result$p_value,
        FDR = NA, Significance = NA, stringsAsFactors = FALSE
      ))
    }
  }

  results$FDR <- p.adjust(results$PValue, method = "BH")
  results$Significance <- ifelse(results$FDR < fdr_threshold, "Significant", "NS")
  results[order(results$FDR), ]
}

Fisher_results <- perform_fisher_tests(DEGs_list, TF_list, Background)

# ---- Restrict to TFs that are themselves significantly differentially active (FDR < 0.05), merge log2FC ----
Fisher_results$DF_TF <- Fisher_results$TF_symbol %in% sig_TFs
Fisher_Merged_TF_DF <- merge(Fisher_results, DEG_TFs, by = "TF_symbol")

max_Odd_plus_ONE <- max(Fisher_Merged_TF_DF$OddsRatio[is.finite(Fisher_Merged_TF_DF$OddsRatio)]) + 1
Fisher_Merged_TF_DF$OddsRatio[is.infinite(Fisher_Merged_TF_DF$OddsRatio)] <- max_Odd_plus_ONE

Fisher_Disreg <- Fisher_Merged_TF_DF %>%
  filter(DF_TF, Significance == "Significant") %>%
  arrange(log2FC)
Fisher_Disreg$Set <- factor(Fisher_Disreg$TF_name, levels = unique(Fisher_Disreg$TF_name))

# ---- Panel G: log2FC (y) x TF (x), size = OddsRatio, color = log2FC ----
p_G_dotplot <- ggplot(Fisher_Disreg, aes(x = Set, y = log2FC, size = OddsRatio, color = log2FC)) +
  geom_point() +
  theme_bw() +
  RotatedAxis() +
  scale_color_gradientn(
    colors = c("blue", "lightgrey", "red"),
    values = scales::rescale(c(min(Fisher_Disreg$log2FC), 0, max(Fisher_Disreg$log2FC)), midpoint = 0)
  ) +
  scale_size_continuous(range = c(0.1, 2)) +
  ylim(-0.2, 0.2) +
  theme(
    text = element_text(family = "Helvetica", size = 6),
    legend.position = "right",
    legend.box = "horizontal",
    axis.text.x = element_text(angle = 90, hjust = 1)
  ) +
  guides(color = guide_colorbar(title.position = "top"), size = guide_legend(title.position = "top"))

save_table_safe(Fisher_Merged_TF_DF, "Figure2_G_Fisher_TF_enrichment_ExtGenes_full.csv")
save_table_safe(Fisher_Disreg,       "Figure2_G_Fisher_TF_enrichment_ExtGenes_significant.csv")

width_in  <- 170 * 0.03937
height_in <- 32  * 0.03937
save_plot_safe(p_G_dotplot, "Figure2_G_TF_regulon_dotplot_ExtGenes", width = width_in, height = height_in)

cat("Significant, differentially-active TFs plotted:", nrow(Fisher_Disreg), "\n")
