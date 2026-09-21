# ============================================================================
# Figure1_00_build_pseudobulk.R
# Project: Decoding Immune Landscapes and Therapeutic Targets in Intracerebral Hemorrhage
# Author:  Dimitrios Kyriakis
# ============================================================================

# ==============================================================================
# Figure 1 — pseudobulk aggregation (feeds panels C/D/E/F)
#
# Loads the Immune Seurat object, aggregates to pseudobulk by CellType x donor,
# and runs dreamlet::processAssays() for normalization. This is the shared
# input every other Figure1_*.R script depends on, so it is computed once here
# and cached to DERIVED_DIR.
#
# Input:  2024_11_Final_Immune.rds
# Output: DERIVED_DIR/Figure1_pbObj_Immune.rds
#         DERIVED_DIR/Figure1_res_proc.rds
#         FIG_DIR/qc_voom_diagnostics.png
#         FIG_DIR/qc_variance_partition.png
# ==============================================================================

source("../../00_setup.R")

Immune <- readRDS(file.path(BASE_DIR, "2024_11_Final_Immune.rds"))

# ---- Convert to SingleCellExperiment ----
mat <- as.matrix(Immune@assays$RNA@counts)
sce <- SingleCellExperiment(
  assays  = list(counts = mat),
  colData = Immune@meta.data
)

# ---- Aggregate to pseudobulk (CellType x donor) ----
pbObj_Immune <- aggregateToPseudoBulk(
  sce,
  assay      = "counts",
  cluster_id = "CellType",
  sample_id  = "donor",
  BPPARAM    = SnowParam(5, progressbar = TRUE)
)

# ---- Normalize / voom-transform each assay ----
res.proc <- processAssays(
  pbObj_Immune,
  ~ (1 | sex) + (1 | dx) + (1 | BL) + age + TSH + mRS_binned,
  min.count = 5,
  BPPARAM   = SnowParam(20, progressbar = TRUE)
)

saveRDS(pbObj_Immune, file.path(DERIVED_DIR, "Figure1_pbObj_Immune.rds"))
saveRDS(res.proc,     file.path(DERIVED_DIR, "Figure1_res_proc.rds"))
cat("Saved: Figure1_pbObj_Immune.rds, Figure1_res_proc.rds\n")

# ------------------------------------------------------------------------------
# QC: voom mean-variance diagnostic + variance partition (not a lettered panel,
# kept for methods/supplementary sanity-checking of the normalization step)
# ------------------------------------------------------------------------------
source("../../00_helpers_stats.R")

p_voom <- my_VroomPlot(res.proc, ncol = 5) + theme(text = element_text(size = 30))
save_plot_safe(p_voom, "qc_voom_diagnostics", width = 16, height = 8)

res.proc_subset <- res.proc[c(2, 3, 4, 5)]
vp.lst <- fitVarPart(res.proc_subset, ~ race + dx + sex + BL + TSH + age + mRS_binned)
p_varpart <- plotVarPart(vp.lst, label.angle = 60, ncol = 5)
save_plot_safe(p_varpart, "qc_variance_partition", width = 16, height = 6)
