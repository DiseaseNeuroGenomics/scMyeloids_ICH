# ============================================================================
# Figure2_analysis_QC_crosscheck.R
# Project: Decoding Immune Landscapes and Therapeutic Targets in Intracerebral Hemorrhage
# Author:  Dimitrios Kyriakis
# ============================================================================

# ==============================================================================
# Figure 2 — cross-figure analysis (not a figure panel): are the MTC_456
# marker genes enriched in the Myeloid lipid/lipoprotein-biosynthesis GO:BP
# terms that were *suppressed* (favorable-outcome) in Figure 1's population
# -level Myeloid GSEA? Hypergeometric overlap test.
#
# Reproduces the notebook's own final resolution (cell 81/82): DEGs1
# (FindAllMarkers, logfc.threshold=0) is the correct gene source per the
# Methods, not the exploratory deg_pa (logfc.threshold=0.25) pass from
# Figure2_E_compareCluster_enrichment.R.
#
# Input:  DERIVED_DIR/Figure1_leading_edge_list.rds        (Figure1_analysis_QC_crosscheck.R)
#         TABLE_DIR/Figure2_master_DEGs_FindAllMarkers.csv (Figure2_master_DEGs.R)
#         TABLE_DIR/Figure1_master_DEGs_all_assays.csv     (Figure1_master_DEGs.R)
# Output: printed hypergeometric test results (console only)
# ==============================================================================

source("../../00_setup.R")

leading_edge_list <- readRDS(file.path(DERIVED_DIR, "Figure1_leading_edge_list.rds"))
DEGs1     <- read.csv(file.path(TABLE_DIR, "Figure2_master_DEGs_FindAllMarkers.csv"))
DEGs_file <- read.csv(file.path(TABLE_DIR, "Figure1_master_DEGs_all_assays.csv"))

cat("MTC_456 DEGs (p_val_adj<0.01):",
    DEGs1 %>% filter(p_val_adj < 0.01, cluster == "MTC_456") %>% nrow(), "\n")
cat("MTC_123 DEGs (p_val_adj<0.01):",
    DEGs1 %>% filter(p_val_adj < 0.01, cluster == "MTC_123") %>% nrow(), "\n")

mtc456_genes <- DEGs1 %>%
  filter(cluster == "MTC_456", p_val_adj < 0.01) %>%
  pull(gene)

overlap_check <- lapply(leading_edge_list, function(genes) intersect(genes, mtc456_genes))
print(overlap_check)

total_genes_tested <- length(unique(DEGs_file$ID[DEGs_file$assay == "Myeloid"]))

for (term in names(leading_edge_list)) {
  le_genes <- leading_edge_list[[term]]
  overlap_n <- length(intersect(le_genes, mtc456_genes))
  p_val <- phyper(
    q = overlap_n - 1,
    m = length(mtc456_genes),
    n = total_genes_tested - length(mtc456_genes),
    k = length(le_genes),
    lower.tail = FALSE
  )
  cat(sprintf("%s: %d/%d overlap, hypergeometric p = %.4f\n", term, overlap_n, length(le_genes), p_val))
}

# ---- Reproducibility record ----
writeLines(capture.output(sessionInfo()), file.path(TABLE_DIR, "Figure2_sessionInfo.txt"))
cat("Saved: Figure2_sessionInfo.txt\n")
