# ============================================================================
# Figure2_D_heatmap_top_genes.R
# Project: Decoding Immune Landscapes and Therapeutic Targets in Intracerebral Hemorrhage
# Author:  Dimitrios Kyriakis
# ============================================================================

# ==============================================================================
# Figure 2, Panel D — z-scored marker-gene heatmap across the 7 MTC subclasses,
# top DE genes per MTC_123/MTC_456 group with key genes labeled
#
# Input:  DERIVED_DIR/Figure2_Metacells.rds                          (Figure2_A_umap_metacells.R)
#         TABLE_DIR/Figure2_master_DEGs_FindAllMarkers.csv           (Figure2_master_DEGs.R)
# Output: FIG_DIR/Figure2_D_heatmap.(pdf|png)
# ==============================================================================

source("../../setup/00_setup.R")
source("../../setup/00_theme_colors.R")

Metacells <- readRDS(file.path(DERIVED_DIR, "Figure2_Metacells.rds"))
DEGs1 <- read.csv(file.path(TABLE_DIR, "Figure2_master_DEGs_FindAllMarkers.csv"))

subset_two_cl <- subset(Metacells, subset = TwoClusters %in% c("MTC_123", "MTC_456"))

top15 <- DEGs1 %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 0.2, p_val_adj < 0.05) %>%
  arrange(desc(avg_log2FC)) %>%
  ungroup()

AverageExpression_RNA <- AverageExpression(
  subset_two_cl, group.by = c("Clusters", "TwoClusters"), features = top15$gene, return.seurat = TRUE
)

label <- as.vector(colnames(AverageExpression_RNA))
split_strings <- strsplit(label, "_")
before_underscore <- sapply(split_strings, function(x) x[1])
after_underscore  <- sapply(split_strings, function(x) x[2])

# AverageExpression(group.by = c("Clusters","TwoClusters")) sanitizes any "_"
# already inside the factor levels to "-" before joining columns with its own
# "_" separator (e.g. "MTC_1" + "MTC_123" -> "MTC-1_MTC-123"), so
# before_underscore comes back as "MTC-1" etc. Un-sanitize it to match
# MTC_Colors' underscore naming (after_underscore/"MTC-123" already matches
# TwoClass_Colors as-is, left untouched).
before_underscore <- gsub("-", "_", before_underscore)

gene_expression_matrix <- LayerData(AverageExpression_RNA, assay = "RNA", layer = "scale.data")
colnames(gene_expression_matrix) <- before_underscore
selected_genes_expression <- as.matrix(gene_expression_matrix[top15$gene, , drop = FALSE])

selected_genes_plot_list <- c(
  "IL3RA", "SPP1", "KLF2", "FOSL2", "FOSL1", "CD163", "MAFB",
  "PPARG", "TCF4", "BHLHE41", "CX3CR1", "IL1B",
  "PTPRG", "ITGA4", "ITGB1", "LILRB1", "LGALS1", "CD63", "CD81", "CD9",
  "SORL1", "APOE", "TREM2"
)
selected_genes_plot <- selected_genes_plot_list[selected_genes_plot_list %in% rownames(selected_genes_expression)]
index_genes <- rownames(selected_genes_expression) %in% selected_genes_plot
selected_genes_plot_ord <- rownames(selected_genes_expression)[index_genes]

top_annotation <- HeatmapAnnotation(
  Clusters  = before_underscore,
  TClusters = after_underscore,
  col = list(TClusters = TwoClass_Colors, Clusters = MTC_Colors)
)

ht <- Heatmap(
  selected_genes_expression, name = "z-score",
  cluster_columns = TRUE, row_km = 3,
  col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
  top_annotation = top_annotation, use_raster = FALSE, show_row_names = TRUE
) +
  rowAnnotation(link = anno_mark(
    at = which(index_genes), labels = selected_genes_plot_ord,
    labels_gp = gpar(fontsize = 11), padding = unit(1, "mm")
  ))

pdf(file.path(FIG_DIR, "Figure2_D_heatmap.pdf"), width = 9, height = 9)
draw(ht)
dev.off()
png(file.path(FIG_DIR, "Figure2_D_heatmap.png"), width = 9, height = 9, units = "in", res = 600)
draw(ht)
dev.off()
cat("Saved: Figure2_D_heatmap (pdf + png)\n")
