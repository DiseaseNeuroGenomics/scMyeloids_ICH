# ============================================================================
# Figure2_supp_enrichment_gsea_gobp_gomf.R
# Project: Decoding Immune Landscapes and Therapeutic Targets in Intracerebral Hemorrhage
# Author:  Dimitrios Kyriakis
# ============================================================================

# ==============================================================================
# Figure 2 (supplementary, not a lettered main-figure panel) — ranked GSEA of
# MTC_456 vs MTC_123 (GO:BP, GO:MF), ranked by the pooled-variance corrected
# t-statistic, with companion ORA. The Hallmark/mTORC1 result is main-text
# Figure 2H — see Figure2_H_hallmark_gsea_mTORC1.R.
#
# Input:  TABLE_DIR/Figure2_master_clean_DEGs_MTC456_vs_MTC123.csv  (Figure2_master_DEGs.R)
# Output: TABLE_DIR/Figure2_supp_ORA_GOBP.csv, _GSEA_GOBP.csv
#         TABLE_DIR/Figure2_supp_ORA_GOMF.csv, _GSEA_GOMF.csv
#         FIG_DIR/Figure2_supp_{ORA,GSEA}_{GOBP,GOMF}.(pdf|png)
# ==============================================================================

source("../../setup/00_setup.R")
source("../../setup/00_helpers_stats.R")

clean_DEGs <- read.csv(file.path(TABLE_DIR, "Figure2_master_clean_DEGs_MTC456_vs_MTC123.csv"))

lfc_cutoff <- 0.5
n_top <- 20

mtc_all <- clean_DEGs %>%
  filter(!is.na(gene), gene != "", !is.na(tstat_corrected), !is.na(log2FC))

sig_genes <- mtc_all %>% filter(p_val_adj < 0.05, abs(log2FC) > lfc_cutoff) %>% pull(gene) %>% unique()
background_genes <- unique(mtc_all$gene)

mtc_ranked <- mtc_all %>%
  group_by(gene) %>%
  summarise(tstat_corrected = mean(tstat_corrected), .groups = "drop") %>%
  arrange(desc(tstat_corrected))
geneList <- sort(setNames(mtc_ranked$tstat_corrected, mtc_ranked$gene), decreasing = TRUE)

top_sig <- function(df, n = n_top) {
  df_sig <- df %>% filter(p.adjust < 0.05)
  if (nrow(df_sig) < n) cat("Note:", nrow(df_sig), "terms pass FDR<0.05 - fewer than", n, "requested.\n")
  df_sig %>% arrange(p.adjust) %>% head(n)
}

run_ora_gsea <- function(gmt_filename, prefix_strip, label) {
  term2gene <- read.gmt(file.path(GMT_DIR, gmt_filename))
  clean_term <- function(x) gsub("_", " ", tolower(gsub(paste0("^", prefix_strip, "_"), "", x)))

  ora <- enricher(gene = sig_genes, universe = background_genes, TERM2GENE = term2gene,
                   pAdjustMethod = "BH", pvalueCutoff = 1, qvalueCutoff = 1)
  ora_df <- as.data.frame(ora)
  ora_df$GeneRatioNum <- sapply(ora_df$GeneRatio, function(x) eval(parse(text = x)))

  set.seed(12345)
  gsea <- GSEA(geneList = geneList, TERM2GENE = term2gene, minGSSize = 15, maxGSSize = 500,
               pvalueCutoff = 1, eps = 0, seed = TRUE, nPermSimple = 10000)
  gsea_df <- as.data.frame(gsea)

  ora_top <- top_sig(ora_df) %>%
    mutate(Description = clean_term(Description),
           Description = factor(Description, levels = Description[order(GeneRatioNum)]))
  p_ora <- ggplot(ora_top, aes(x = GeneRatioNum, y = Description)) +
    geom_point(aes(size = Count, color = p.adjust)) +
    scale_color_gradient(low = "#B2182B", high = "#4C78A8") +
    labs(x = "Gene Ratio", y = NULL, title = paste0("MTC_456 vs MTC_123 ORA - ", label, ", top ", nrow(ora_top), " (FDR<0.05)"),
         color = "FDR", size = "Gene Count") +
    theme_classic(base_size = 11) + theme(axis.text = element_text(color = "black")) + theme_natmed()

  gsea_top <- top_sig(gsea_df) %>%
    mutate(Description = clean_term(Description),
           Description = factor(Description, levels = Description[order(NES)]))
  p_gsea <- ggplot(gsea_top, aes(x = NES, y = Description)) +
    geom_point(aes(size = -log10(p.adjust), color = NES)) +
    scale_color_gradient2(low = "#4C78A8", mid = "white", high = "#D95F5F", midpoint = 0) +
    geom_vline(xintercept = 0, linewidth = 0.4, linetype = "dashed", color = "grey40") +
    labs(x = "Normalized Enrichment Score", y = NULL,
         title = paste0("MTC_456 vs MTC_123 GSEA - ", label, ", top ", nrow(gsea_top), " (FDR<0.05)"),
         color = "NES", size = expression(-log[10] ~ (FDR))) +
    theme_classic(base_size = 11) + theme(axis.text = element_text(color = "black")) + theme_natmed()

  list(ora_df = ora_df, gsea_df = gsea_df, p_ora = p_ora, p_gsea = p_gsea)
}

res_gobp <- run_ora_gsea("c5.go.bp.v7.5.1.symbols.gmt", "GOBP", "GO:BP")
save_table_safe(res_gobp$ora_df,  "Figure2_supp_ORA_GOBP.csv")
save_table_safe(res_gobp$gsea_df, "Figure2_supp_GSEA_GOBP.csv")
save_plot_safe(res_gobp$p_ora,  "Figure2_supp_ORA_GOBP",  width = 7, height = 5)
save_plot_safe(res_gobp$p_gsea, "Figure2_supp_GSEA_GOBP", width = 7, height = 5)

res_gomf <- run_ora_gsea("c5.go.mf.v7.5.1.symbols.gmt", "GOMF", "GO:MF")
save_table_safe(res_gomf$ora_df,  "Figure2_supp_ORA_GOMF.csv")
save_table_safe(res_gomf$gsea_df, "Figure2_supp_GSEA_GOMF.csv")
save_plot_safe(res_gomf$p_ora,  "Figure2_supp_ORA_GOMF",  width = 7, height = 5)
save_plot_safe(res_gomf$p_gsea, "Figure2_supp_GSEA_GOMF", width = 7, height = 5)
