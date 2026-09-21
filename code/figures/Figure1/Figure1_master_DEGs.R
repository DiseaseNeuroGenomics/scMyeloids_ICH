# ============================================================================
# Figure1_master_DEGs.R
# Project: Decoding Immune Landscapes and Therapeutic Targets in Intracerebral Hemorrhage
# Author:  Dimitrios Kyriakis
# ============================================================================

# ==============================================================================
# Figure 1 — master differential expression table (feeds panels D/E/F)
#
# dreamlet DEG test (mRS poor vs good outcome) across all cell types. This is
# the single master table that Figure1_E_deg_summary_bar.R and
# Figure1_F_myeloid_enrichment.R both load, rather than each recomputing it.
#
# Input:  DERIVED_DIR/Figure1_res_proc.rds  (from Figure1_00_build_pseudobulk.R)
# Output: TABLE_DIR/Figure1_master_DEGs_all_assays.csv
# ==============================================================================

source("../../00_setup.R")

res.proc <- readRDS(file.path(DERIVED_DIR, "Figure1_res_proc.rds"))

res.dream_mRS <- dreamlet(
  res.proc, ~ (1 | race) + (1 | dx) + (1 | BL) + age + mRS_binned + 0,
  contrasts = c(Diff_B_vs_G = "mRS_binnedB-mRS_binnedG"),
  BPPARAM = SnowParam(20, progressbar = TRUE)
)

DEGs_file <- as.data.frame(topTable(
  res.dream_mRS, p.value = 1, lfc = 0, number = 100000000, coef = "Diff_B_vs_G"
))

DEGs_file <- DEGs_file[order(DEGs_file$adj.P.Val), ]
DEGs_file$cluster <- ifelse(DEGs_file$logFC > 0, "Unfavorable", "Favorable")
DEGs_file$cluster <- paste0(DEGs_file$cluster, "_", DEGs_file$assay)
DEGs_file$gene <- DEGs_file$ID
DEGs_file <- DEGs_file %>% group_by(ID) %>% mutate(pos_avgLFC = mean(logFC)) %>% ungroup()

save_table_safe(DEGs_file, "Figure1_master_DEGs_all_assays.csv")
