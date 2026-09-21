# ============================================================================
# Figure2_E_compareCluster_enrichment.R
# Project: Decoding Immune Landscapes and Therapeutic Targets in Intracerebral Hemorrhage
# Author:  Dimitrios Kyriakis
# ============================================================================

# ==============================================================================
# Figure 2, Panel E — GO (MF + BP) over-representation dotplot, MTC_123 vs
# MTC_456 side by side (clusterProfiler::compareCluster)
#
# Note: uses Seurat's default logfc.threshold (0.25), unlike the
# logfc.threshold=0 master DEG table in Figure2_master_DEGs.R — this is a
# separate significant-gene call specific to this panel's ORA, not a rerun of
# the master table.
#
# Input:  DERIVED_DIR/Figure2_Metacells.rds  (from Figure2_A_umap_metacells.R)
# Output: TABLE_DIR/Figure2_E_deg_pa_FindAllMarkers.csv
#         FIG_DIR/Figure2_E_GO_MF_dotplot.(pdf|png)
#         FIG_DIR/Figure2_E_GO_BP_dotplot.(pdf|png)
# ==============================================================================

source("../../00_setup.R")
library(org.Hs.eg.db)
library(enrichplot)

Metacells <- readRDS(file.path(DERIVED_DIR, "Figure2_Metacells.rds"))

Metacells_2 <- subset(Metacells, subset = TwoClusters %in% c("MTC_123", "MTC_456"))
Idents(Metacells_2) <- "TwoClusters"

deg_pa <- FindAllMarkers(Metacells_2, only.pos = TRUE)
significant_Genes <- deg_pa %>% filter(p_val_adj < 0.01, avg_log2FC > 1.5)

formula_MF <- compareCluster(
  gene ~ cluster, data = significant_Genes, fun = "enrichGO",
  OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont = "MF",
  pAdjustMethod = "BH", pvalueCutoff = 0.05
)
formula_BP <- compareCluster(
  gene ~ cluster, data = significant_Genes, fun = "enrichGO",
  OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont = "BP",
  pAdjustMethod = "BH", pvalueCutoff = 0.05
)

p_E_mf <- dotplot(formula_MF, showCategory = 10) +
  theme_classic(base_size = 14) +
  theme(axis.text.y = element_text(color = "black"))

p_E_bp <- dotplot(formula_BP, showCategory = 10) +
  theme_classic(base_size = 14) +
  theme(axis.text.y = element_text(color = "black"))

save_table_safe(deg_pa, "Figure2_E_deg_pa_FindAllMarkers.csv")
save_table_safe(formula_MF, "Figure2_E_GO_MF_dotplot.csv")
save_table_safe(formula_BP, "Figure2_E_GO_BP_dotplot.csv")
save_plot_safe(p_E_mf, "Figure2_E_GO_MF_dotplot", width = 7, height = 10)
save_plot_safe(p_E_bp, "Figure2_E_GO_BP_dotplot", width = 7, height = 10)
