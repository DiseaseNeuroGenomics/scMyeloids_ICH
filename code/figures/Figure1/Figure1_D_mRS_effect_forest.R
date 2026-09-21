# ============================================================================
# Figure1_D_mRS_effect_forest.R
# Project: Decoding Immune Landscapes and Therapeutic Targets in Intracerebral Hemorrhage
# Author:  Dimitrios Kyriakis
# ============================================================================

# ==============================================================================
# Figure 1, Panel D — mRS effect size forest plot per cell type
#
# crumblr variance partition + dream() mixed-model contrast (mRS poor vs good
# outcome), meta-analyzed per cell type and plotted against a composition
# dendrogram.
#
# Input:  DERIVED_DIR/Figure1_pbObj_Immune.rds  (from Figure1_00_build_pseudobulk.R)
# Output: TABLE_DIR/Figure1_D_res_mRS_topTable.csv
#         TABLE_DIR/Figure1_D_res_meta_mRS.csv
#         FIG_DIR/Figure1_D_mRS_effect_forest.(pdf|png)
# ==============================================================================

source("../../00_setup.R")
source("../../00_helpers_stats.R")

pbObj_Immune <- readRDS(file.path(DERIVED_DIR, "Figure1_pbObj_Immune.rds"))

# ---- CrumblR: cell composition variance partition ----
cobj <- crumblr(cellCounts(pbObj_Immune))
vp <- fitExtractVarPartModel(cobj, ~ (1 | race) + (1 | sex) + (1 | BL) + (1 | dx) + age + TSH,
                              colData(pbObj_Immune))
p_composition_varpart <- plotPercentBars(vp)
save_plot_safe(p_composition_varpart, "qc_crumblr_composition_variance", width = 6, height = 4.5)

# ---- Composition dendrogram ----
hc <- buildClusterTreeFromPB(pbObj_Immune)
fig.tree <- dh_plotTree(ape::as.phylo(hc), xmax.scale = 2.5) +
  theme(legend.position = "right", text = element_text(size = 15))

# ---- mRS mixed-model contrast: mRS_binned "Bad" vs "Good" ----
form <- ~ (1 | race) + (1 | sex) + (1 | BL) + (1 | dx) + age + mRS_binned + 0
L <- makeContrastsDream(form, colData(pbObj_Immune),
                         contrasts = c(mRS_contrast = "mRS_binnedB - mRS_binnedG"))
fit <- dream(cobj, form, colData(pbObj_Immune), L = L)
fit <- eBayes(fit)

res.mRS <- topTable(fit, coef = "mRS_contrast", number = Inf, sort.by = "none")
res.mRS$assay <- rownames(res.mRS)

res.meta_mRS <- dh_meta_analysis(list(res.mRS))

res.meta_mRS$FDR <- p.adjust(res.meta_mRS$p.value, method = "BH")
res.meta_mRS$log10FDR <- -log10(res.meta_mRS$FDR)

p_D_forest <- dh_plotCoef(res.meta_mRS, fig.tree, ylab = "mRS Effect Size (logFC)") +
  theme_natmed() +
  theme(legend.position = c(1.23, 0.25),
        legend.justification = c("right", "top"),
        legend.background = element_rect(fill = alpha("white", 0.6), color = NA),
        legend.key = element_rect(fill = alpha("white", 0)))

save_table_safe(res.mRS,      "Figure1_D_res_mRS_topTable.csv")
save_table_safe(res.meta_mRS, "Figure1_D_res_meta_mRS.csv")
save_plot_safe(p_D_forest, "Figure1_D_mRS_effect_forest", width = 4, height = 4)
