# ============================================================================
# Figure2_H_hallmark_gsea_mTORC1.R
# Project: Decoding Immune Landscapes and Therapeutic Targets in Intracerebral Hemorrhage
# Author:  Dimitrios Kyriakis
# ============================================================================

# ==============================================================================
# Figure 2, Panel H — Hallmark GSEA of MTC_456 vs MTC_123, ranked by the
# pooled-variance corrected t-statistic, with HALLMARK_MTORC1_SIGNALING
# outlined (cited in Discussion as directly linking the MTC_456 transcriptional
# state to mTOR pathway activation).
#
# Input:  TABLE_DIR/Figure2_master_clean_DEGs_MTC456_vs_MTC123.csv  (Figure2_master_DEGs.R)
#         GMT_DIR/h.all.v7.5.1.symbols.gmt
# Output: TABLE_DIR/Figure2_H_GSEA_Hallmark_full.csv
#         FIG_DIR/Figure2_H_hallmark_gsea_mTORC1.(pdf|png)
# ==============================================================================

source("../../00_setup.R")
source("../../00_helpers_stats.R")

clean_DEGs <- read.csv(file.path(TABLE_DIR, "Figure2_master_clean_DEGs_MTC456_vs_MTC123.csv"))

n_top <- 20

mtc_all <- clean_DEGs %>%
  filter(!is.na(gene), gene != "", !is.na(tstat_corrected), !is.na(log2FC))

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

term2gene_hallmark <- read.gmt(file.path(GMT_DIR, "h.all.v7.5.1.symbols.gmt"))
clean_hallmark <- function(x) gsub("_", " ", tolower(gsub("^HALLMARK_", "", x)))

set.seed(12345)
gsea_hallmark <- GSEA(geneList = geneList, TERM2GENE = term2gene_hallmark,
                       minGSSize = 15, maxGSSize = 500, pvalueCutoff = 1,
                       eps = 0, seed = TRUE, nPermSimple = 10000)
gsea_df_hallmark <- as.data.frame(gsea_hallmark) %>% arrange(p.adjust) %>% mutate(rank = row_number())

mtor_row <- gsea_df_hallmark %>% filter(grepl("MTORC1", ID, ignore.case = TRUE))
cat("HALLMARK_MTORC1_SIGNALING result:\n")
print(mtor_row[, c("ID", "NES", "pvalue", "p.adjust", "setSize", "rank")])

gsea_top_hallmark <- top_sig(gsea_df_hallmark) %>%
  mutate(is_mtor = grepl("MTORC1", ID, ignore.case = TRUE),
         Description = clean_hallmark(Description),
         Description = factor(Description, levels = Description[order(NES)]))

p_H_hallmark <- ggplot(gsea_top_hallmark, aes(x = NES, y = Description)) +
  geom_point(aes(size = -log10(p.adjust), fill = NES), shape = 21, color = "grey30", stroke = 0.3) +
  geom_point(data = filter(gsea_top_hallmark, is_mtor),
             aes(size = -log10(p.adjust)), shape = 21, fill = NA, color = "black", stroke = 1.5) +
  scale_fill_gradient2(low = "#4C78A8", mid = "white", high = "#D95F5F", midpoint = 0) +
  geom_vline(xintercept = 0, linewidth = 0.4, linetype = "dashed", color = "grey40") +
  labs(x = "Normalized Enrichment Score", y = NULL,
       title = paste0("MTC_456 vs MTC_123 GSEA - Hallmark, top ", nrow(gsea_top_hallmark), " (FDR<0.05), mTORC1 outlined bold"),
       fill = "NES", size = expression(-log[10] ~ (FDR))) +
  theme_classic(base_size = 11) + theme(axis.text = element_text(color = "black")) + theme_natmed()

save_table_safe(gsea_df_hallmark, "Figure2_H_GSEA_Hallmark_full.csv")
save_plot_safe(p_H_hallmark, "Figure2_H_hallmark_gsea_mTORC1", width = 5.2, height = 5)
