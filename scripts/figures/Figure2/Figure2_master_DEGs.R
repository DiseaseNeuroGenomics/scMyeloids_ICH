# ============================================================================
# Figure2_master_DEGs.R
# Project: Decoding Immune Landscapes and Therapeutic Targets in Intracerebral Hemorrhage
# Author:  Dimitrios Kyriakis
# ============================================================================

# ==============================================================================
# Figure 2 — master differential expression table (feeds panels D/E and the
# IREA cytokine-enrichment supplementary analysis)
#
# FindAllMarkers(MTC_123 vs MTC_456) plus a corrected (pooled-variance) t-stat,
# since Seurat's default avg_log2FC/p_val_adj don't carry group sample sizes
# needed for a proper ranked GSEA.
#
# Input:  DERIVED_DIR/Figure2_Metacells.rds  (from Figure2_A_umap_metacells.R)
# Output: TABLE_DIR/Figure2_master_DEGs_FindAllMarkers.csv         (DEGs1)
#         TABLE_DIR/Figure2_master_clean_DEGs_MTC456_vs_MTC123.csv (clean_DEGs)
# ==============================================================================

source("../../setup/00_setup.R")

Metacells <- readRDS(file.path(DERIVED_DIR, "Figure2_Metacells.rds"))

subset_two_cl <- subset(Metacells, subset = TwoClusters %in% c("MTC_123", "MTC_456"))
Idents(subset_two_cl) <- "TwoClusters"

DEGs1 <- FindAllMarkers(subset_two_cl, only.pos = TRUE, logfc.threshold = 0, return.thresh = 1)

# ---- Corrected (pooled-variance) t-statistic, for ranked GSEA ----
mat_456 <- subset_two_cl@assays$RNA@data[, subset_two_cl$TwoClusters == "MTC_456"]
mat_123 <- subset_two_cl@assays$RNA@data[, subset_two_cl$TwoClusters == "MTC_123"]
n_456 <- ncol(mat_456)
n_123 <- ncol(mat_123)

mean_456 <- Matrix::rowMeans(mat_456)
mean_123 <- Matrix::rowMeans(mat_123)
var_456 <- (Matrix::rowMeans(mat_456^2) - mean_456^2) * (n_456 / (n_456 - 1))
var_123 <- (Matrix::rowMeans(mat_123^2) - mean_123^2) * (n_123 / (n_123 - 1))
se_corrected <- sqrt((var_456 / n_456) + (var_123 / n_123))
tstat_corrected <- (mean_456 - mean_123) / se_corrected

fold_ch_res <- FoldChange(
  subset_two_cl, ident.1 = "MTC_456", ident.2 = "MTC_123",
  logfc.threshold = 0, slot = "data", min.pct = 0.1, min.cells.group = 3,
  pseudocount.use = 1, mean.fxn = NULL, fc.name = NULL, base = 2, return.thresh = 1
)

# sd.1/sd.2 reuse var_456/var_123 computed above instead of apply(mat, 1, sd):
# apply() on a sparse dgCMatrix silently coerces the whole matrix to dense
# first, which is what was blowing past the job's memory allocation. Same
# Bessel-corrected variance, computed sparsity-safe, no extra pass needed.
fold_ch_res$sd.1 <- sqrt(var_456)
fold_ch_res$sd.2 <- sqrt(var_123)

averageMT <- as.data.frame(AverageExpression(subset_two_cl, group.by = "TwoClusters"))

merged_df <- merge(DEGs1, fold_ch_res, by = "row.names", all = TRUE)
rownames(merged_df) <- merged_df$Row.names
merged_df[, 1] <- NULL

merged_df2 <- merge(merged_df, averageMT, by = "row.names", all = TRUE)
rownames(merged_df2) <- merged_df2$Row.names
merged_df2[, 1] <- NULL

merged_df2$se <- sqrt((merged_df2$sd.1^2 / n_456) + (merged_df2$sd.2^2 / n_123))
merged_df2$log2FC <- ifelse(merged_df2$cluster == "MTC_456", merged_df2$avg_log2FC.x, -1 * merged_df2$avg_log2FC.x)
merged_df2$tstat <- merged_df2$log2FC / merged_df2$se
merged_df2$tstat_corrected <- tstat_corrected[rownames(merged_df2)]

clean_DEGs <- merged_df2 %>% na.omit()
clean_DEGs$gene <- rownames(clean_DEGs)
clean_DEGs <- clean_DEGs %>% arrange(desc(log2FC))

cat("cor(|tstat|, |tstat_corrected|) =", cor(abs(clean_DEGs$tstat), abs(clean_DEGs$tstat_corrected)), "\n")

save_table_safe(DEGs1,      "Figure2_master_DEGs_FindAllMarkers.csv")
save_table_safe(clean_DEGs, "Figure2_master_clean_DEGs_MTC456_vs_MTC123.csv")
readr::write_tsv(clean_DEGs, file.path(TABLE_DIR, "Figure2_master_clean_DEGs_MTC456_vs_MTC123_full_incl_tstat.tsv"))
