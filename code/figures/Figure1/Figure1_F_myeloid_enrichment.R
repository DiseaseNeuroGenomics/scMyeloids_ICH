# ============================================================================
# Figure1_F_myeloid_enrichment.R
# Project: Decoding Immune Landscapes and Therapeutic Targets in Intracerebral Hemorrhage
# Author:  Dimitrios Kyriakis
# ============================================================================

# ==============================================================================
# Figure 1, Panel F — Myeloid GO:BP enrichment (GSEA, ranked by log2FC), plus
# the ORA companion analysis
#
# Input:  TABLE_DIR/Figure1_master_DEGs_all_assays.csv  (from Figure1_master_DEGs.R)
#         GMT_DIR/c5.go.bp.v7.5.1.symbols.gmt
# Output: TABLE_DIR/Figure1_F_myeloid_ORA_full.csv, _top20.csv
#         TABLE_DIR/Figure1_F_myeloid_GSEA_full.csv, _top20.csv
#         FIG_DIR/Figure1_F_myeloid_gsea_gobp.(pdf|png)
#         FIG_DIR/Figure1_F_myeloid_ora_gobp.(pdf|png)   (supplementary companion)
# ==============================================================================

source("../../00_setup.R")
source("../../00_helpers_stats.R")

DEGs_file <- read.csv(file.path(TABLE_DIR, "Figure1_master_DEGs_all_assays.csv"))

gmt_file  <- file.path(GMT_DIR, "c5.go.bp.v7.5.1.symbols.gmt")
term2gene <- read.gmt(gmt_file)

lfc_cutoff <- 0.5
n_top <- 20

clean_term <- function(x) tolower(gsub("_", " ", gsub("^GOBP_", "", x)))

myeloid_all <- DEGs_file %>%
  filter(assay == "Myeloid", !is.na(gene), gene != "", !is.na(t), !is.na(logFC))

# ---- ORA: significant genes vs tested-gene background ----
sig_genes <- myeloid_all %>%
  filter(adj.P.Val < 0.05, abs(logFC) > lfc_cutoff) %>%
  pull(gene) %>% unique()
background_genes <- unique(myeloid_all$gene)

ora_myeloid <- enricher(
  gene = sig_genes, universe = background_genes, TERM2GENE = term2gene,
  pAdjustMethod = "BH", pvalueCutoff = 1, qvalueCutoff = 1
)
ora_df <- as.data.frame(ora_myeloid)
ora_df$GeneRatioNum <- sapply(ora_df$GeneRatio, function(x) eval(parse(text = x)))

# ---- GSEA: full ranked list, moderated t-statistic ----
myeloid_ranked <- myeloid_all %>%
  group_by(gene) %>%
  summarise(t = mean(t), .groups = "drop") %>%
  arrange(desc(t))
geneList <- sort(setNames(myeloid_ranked$t, myeloid_ranked$gene), decreasing = TRUE)

set.seed(12345)
gsea_myeloid <- GSEA(
  geneList = geneList, TERM2GENE = term2gene,
  minGSSize = 15, maxGSSize = 500, pvalueCutoff = 1, eps = 0, seed = TRUE
)
gsea_df <- as.data.frame(gsea_myeloid)

# ---- ORA plot: top 20 terms by FDR ----
ora_top <- ora_df %>%
  arrange(p.adjust) %>%
  head(n_top) %>%
  mutate(Description = clean_term(Description),
         Description = factor(Description, levels = Description[order(GeneRatioNum)]))

p_F_ora <- ggplot(ora_top, aes(x = GeneRatioNum, y = Description)) +
  geom_point(aes(size = Count, color = p.adjust)) +
  scale_color_gradient(low = "#B2182B", high = "#4C78A8") +
  labs(x = "Gene Ratio", y = NULL, title = paste0("Myeloid ORA - top ", n_top, " terms"),
       color = "FDR", size = "Gene Count") +
  theme_classic(base_size = 11) + theme(axis.text = element_text(color = "black")) +
  theme_natmed()

# ---- GSEA plot (panel F): top 10 up + top 10 down by FDR ----
gsea_top <- bind_rows(
  gsea_df %>% filter(NES > 0) %>% arrange(p.adjust) %>% head(n_top / 2),
  gsea_df %>% filter(NES < 0) %>% arrange(p.adjust) %>% head(n_top / 2)
) %>%
  mutate(Description = clean_term(Description),
         Description = factor(Description, levels = Description[order(NES)]))

p_F_gsea <- ggplot(gsea_top, aes(x = NES, y = Description)) +
  geom_point(aes(size = -log10(p.adjust), color = NES)) +
  scale_color_gradient2(low = "#4C78A8", mid = "white", high = "#D95F5F", midpoint = 0) +
  geom_vline(xintercept = 0, linewidth = 0.4, linetype = "dashed", color = "grey40") +
  labs(x = "Normalized Enrichment Score", y = NULL,
       title = paste0("Myeloid GSEA - top ", n_top, " terms"),
       color = "NES", size = expression(-log[10] ~ (FDR))) +
  theme_classic(base_size = 11) + theme(axis.text = element_text(color = "black")) +
  theme_natmed() +
  theme(panel.background = element_rect(fill = "white", color = NA),
        plot.background  = element_rect(fill = "white", color = NA),
        legend.key       = element_rect(fill = "white", color = NA),
        legend.background = element_rect(fill = "white", color = NA),
        axis.line = element_line(color = "black", linewidth = 0.4))

save_table_safe(ora_df,   "Figure1_F_myeloid_ORA_full.csv")
save_table_safe(gsea_df,  "Figure1_F_myeloid_GSEA_full.csv")
save_table_safe(ora_top,  "Figure1_F_myeloid_ORA_top20.csv")
save_table_safe(gsea_top, "Figure1_F_myeloid_GSEA_top20.csv")

save_plot_safe(p_F_gsea, "Figure1_F_myeloid_gsea_gobp", width = 4, height = 4)
save_plot_safe(p_F_ora,  "Figure1_F_myeloid_ora_gobp",  width = 8, height = 5)
