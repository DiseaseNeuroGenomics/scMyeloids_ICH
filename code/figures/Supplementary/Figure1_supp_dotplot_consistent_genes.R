# ============================================================================
# Figure1_supp_dotplot_consistent_genes.R
# Project: Decoding Immune Landscapes and Therapeutic Targets in Intracerebral Hemorrhage
# Author:  Dimitrios Kyriakis
# ============================================================================

# ==============================================================================
# Figure 1 (supplementary, not a lettered main-figure panel) — dotplot of genes
# significant and directionally consistent across >= N cell-type assays
# (B Cells excluded due to low power).
#
# Input:  TABLE_DIR/Figure1_master_DEGs_all_assays.csv  (from Figure1_master_DEGs.R)
# Output: TABLE_DIR/Figure1_supp_consistent_genes_diagnostic.csv
#         TABLE_DIR/Figure1_supp_consistent_genes_dotplot_data.csv
#         FIG_DIR/Figure1_supp_dotplot_consistent_genes.(pdf|png)
# ==============================================================================

source("../../00_setup.R")
source("../../00_helpers_stats.R")

DEGs_file <- read.csv(file.path(TABLE_DIR, "Figure1_master_DEGs_all_assays.csv"))

# ------------------------------------------------------------------------------
# EDIT: minimum number of assays a gene must be significant+consistent in.
# Re-run after checking the diagnostic counts printed below.
# ------------------------------------------------------------------------------
min_assays_sig <- 2

DEGs_file_noB <- DEGs_file %>% filter(assay != "B Cells")
n_assays <- length(unique(DEGs_file_noB$assay))

gene_summary <- DEGs_file_noB %>%
  filter(adj.P.Val < 0.01, abs(logFC) > 1) %>%
  group_by(ID) %>%
  summarise(n_sig_assays = n_distinct(assay),
            n_pos = sum(logFC > 0),
            n_neg = sum(logFC < 0),
            consistent_direction = (n_pos == n_sig_assays | n_neg == n_sig_assays),
            .groups = "drop")

for (k in 1:n_assays) {
  n_genes <- gene_summary %>% filter(n_sig_assays >= k, consistent_direction) %>% nrow()
  cat("Significant in >=", k, "of", n_assays, "assays (consistent direction):", n_genes, "genes\n")
}

consistent_genes <- gene_summary %>%
  filter(n_sig_assays >= min_assays_sig, consistent_direction) %>%
  pull(ID)
cat("\nUsing min_assays_sig =", min_assays_sig, "->", length(consistent_genes), "genes\n")

sign_DEGS_df <- DEGs_file_noB[DEGs_file_noB$ID %in% consistent_genes, ]

data.plot <- DEGs_file_noB[DEGs_file_noB$ID %in% unique(as.factor(sign_DEGS_df$ID)), ] %>%
  na.omit()
data.plot <- data.plot[order(data.plot$pos_avgLFC), ]
data.plot$gene <- factor(as.vector(data.plot$ID), levels = unique(as.vector(data.plot$ID)))
data.plot$assay <- factor(data.plot$assay, levels = rev(sort(unique(data.plot$assay))))

cols <- c("blue", "red")
p_consistent_dotplot <- ggplot(data.plot, aes(x = gene, y = assay, size = abs(logFC), fill = logFC)) +
  geom_point(shape = 21, colour = "black") +
  guides(size = guide_legend(title = "log2FC"), fill = guide_legend(title = "log2FC")) +
  scale_fill_gradient2(low = cols[1], mid = "white", high = cols[2]) +
  theme_natmed() +
  theme(
    axis.title = element_blank(),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background  = element_rect(fill = "white", color = NA),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
    axis.line = element_line(color = "black", linewidth = 0.4),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    legend.key = element_rect(fill = "white", color = NA),
    legend.background = element_rect(fill = "white", color = NA),
    legend.position = "bottom"
  )

save_table_safe(gene_summary, "Figure1_supp_consistent_genes_diagnostic.csv")
save_table_safe(data.plot,    "Figure1_supp_consistent_genes_dotplot_data.csv")
save_plot_safe(p_consistent_dotplot, "Figure1_supp_dotplot_consistent_genes", width = 13, height = 3)
