# ============================================================================
# Figure2_C_mRS_effect_tree.R
# Project: Decoding Immune Landscapes and Therapeutic Targets in Intracerebral Hemorrhage
# Author:  Dimitrios Kyriakis
# ============================================================================

# ==============================================================================
# Figure 2, Panel C — MTC subclass composition dendrogram with per-node mRS
# significance (treeTest) + mRS effect-size forest plot
#
# Input:  DERIVED_DIR/Figure2_Metacells.rds  (from Figure2_A_umap_metacells.R)
# Output: TABLE_DIR/Figure2_C_res_mRS_topTable.csv
#         TABLE_DIR/Figure2_C_res_meta_mRS.csv
#         FIG_DIR/Figure2_C_tree_treeTest.(pdf|png)
#         FIG_DIR/Figure2_C_coef_forest.(pdf|png)
#         FIG_DIR/Figure2_C_composite_tree_coef.(pdf|png)
# ==============================================================================

source("../../00_setup.R")
source("../../00_helpers_stats.R")

Metacells <- readRDS(file.path(DERIVED_DIR, "Figure2_Metacells.rds"))

# ---- Pseudobulk by MTC subclass x donor ----
mat <- as.matrix(Metacells@assays$RNA@counts)
sce <- SingleCellExperiment(assays = list(counts = mat), colData = Metacells@meta.data)
pbObj_subtype <- aggregateToPseudoBulk(
  sce, assay = "counts", cluster_id = "Clusters", sample_id = "Donor",
  BPPARAM = SnowParam(5, progressbar = TRUE)
)

cobj <- crumblr(cellCounts(pbObj_subtype))
vp <- fitExtractVarPartModel(cobj, ~ (1 | race) + (1 | sex) + (1 | CL) + (1 | dx) + age + TSH,
                              colData(pbObj_subtype))
p_composition_varpart <- plotPercentBars(vp)
save_plot_safe(p_composition_varpart, "qc_crumblr_MTC_composition_variance", width = 6, height = 4.5)

# ---- mRS mixed-model contrast ----
form <- ~ mRS_binned + (1 | sex) + dx + age + CL + 0
L <- makeContrastsDream(form, colData(pbObj_subtype),
                         contrasts = c(mRS_contrast = "mRS_binnedB - mRS_binnedG"))
fit <- dream(cobj, form, colData(pbObj_subtype), L = L)
fit <- eBayes(fit)

res.mRS <- topTable(fit, coef = "mRS_contrast", number = Inf, sort.by = "none")
res.mRS$assay <- rownames(res.mRS)
res.meta_mRS <- dh_meta_analysis(list(res.mRS))

res.meta_mRS$FDR <- p.adjust(res.meta_mRS$p.value, method = "BH")
res.meta_mRS$log10FDR <- -log10(res.meta_mRS$FDR)

hc <- buildClusterTreeFromPB(pbObj_subtype)
fig.tree <- dh_plotTree(ape::as.phylo(hc), xmax.scale = 2.5) +
  theme(legend.position = "bottom", text = element_text(size = 15))

fig.es1 <- dh_plotCoef(res.meta_mRS, fig.tree, ylab = "mRS", label_axis = FALSE)
p_C_composite <- fig.es1 %>% aplot::insert_left(fig.tree, width = 2)

# ---- Per-node significance tree (treeTest) ----
low <- "grey90"; mid <- "red"; high <- "darkred"; xmax.scale <- 3.5

tree_res <- treeTest(fit, cobj, hc, coef = "mRS_contrast")
fig_base <- ggtree(tree_res, branch.length = "none") +
  geom_tiplab(color = "black", size = 4, hjust = 0, offset = 0.2) +
  theme(legend.position = "top left", plot.title = element_text(hjust = 0.5))
xmax <- layer_scales(fig_base)$x$range$range[2]

p_C_tree <- ggtree(tree_res, branch.length = "none") +
  geom_tiplab(color = "black", size = 6, hjust = 0, offset = 0.5) +
  geom_point2(aes(label = node, color = pmin(4, -log10(pvalue)), size = pmin(5, -log10(pvalue)))) +
  scale_color_gradient2(name = bquote(-log[10]~pvalue), limits = c(0, 4), low = low, mid = mid, high = high,
                         midpoint = -log10(0.01)) +
  scale_size_area(name = bquote(-log[10]~pvalue), limits = c(0, 4)) +
  geom_text2(aes(label = '+', subset = pvalue < 0.05), color = "white", size = 6, vjust = 0.3, hjust = 0.5) +
  theme(legend.position = "right", plot.title = element_text(hjust = 0.5)) +
  xlim(0, xmax * xmax.scale) +
  ggtitle("mRS_contrast \n mRS_binnedB -mRS_binnedG")

save_table_safe(res.mRS,      "Figure2_C_res_mRS_topTable.csv")
save_table_safe(res.meta_mRS, "Figure2_C_res_meta_mRS.csv")

save_plot_safe(p_C_tree,       "Figure2_C_tree_treeTest",         width = 6, height = 4.5)
save_plot_safe(fig.es1,        "Figure2_C_coef_forest",           width = 4, height = 4.5)
save_plot_safe(p_C_composite,  "Figure2_C_composite_tree_coef",   width = 6, height = 4.5)
