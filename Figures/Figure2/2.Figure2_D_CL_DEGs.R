#!/usr/bin/env Rscript
# Author: Dimitrios Kyriakis
# Description: Differential gene expression between MTC_456 (poor outcome) and
#              MTC_123 (good outcome) metaclusters using FindAllMarkers.
#              Computes full DEG table with SE/t-stat and writes output CSV/TSV.
#              Produces Figure 2D (ComplexHeatmap of top DEGs).
# Usage: Rscript 2.Figure2_D_CL_DEGs.R <metacell_rds> <output_dir>
#   <metacell_rds>  Path to 2024_03_28_Myeloid_Metacells_Subclass_ADAM.rds.ztsd
#   <output_dir>    Directory for output files (default: current dir)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
    stop("Usage: Rscript 2.Figure2_D_CL_DEGs.R <metacell_rds> [output_dir]")
}

input_rds  <- args[1]
output_dir <- if (length(args) >= 2) args[2] else "."

set.seed(12345)

suppressPackageStartupMessages({
    library(Seurat)
    library(dplyr)
    library(ggplot2)
    library(cowplot)
    library(ComplexHeatmap)
    library(circlize)
    library(RColorBrewer)
})

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

Metacells <- readCRDS(input_rds)

Metacells$TwoClusters <- as.vector(Metacells$Clusters)
Metacells$TwoClusters[Metacells$Clusters %in% c("MTC_1", "MTC_2", "MTC_3")] <- "MTC_123"
Metacells$TwoClusters[Metacells$Clusters %in% c("MTC_4", "MTC_5", "MTC_6")] <- "MTC_456"

Metacells_2 <- subset(Metacells, subset = TwoClusters %in% c("MTC_123", "MTC_456"))
Idents(Metacells_2) <- "TwoClusters"

DEGs <- FindAllMarkers(Metacells_2, only.pos = FALSE)

# Add SE and t-stat
DEGs <- DEGs %>%
    mutate(
        se    = avg_log2FC / ifelse(abs(log2(p_val + 1e-300)), sqrt(log2(p_val + 1e-300) * -2), NA),
        tstat = avg_log2FC / (sd(avg_log2FC) / sqrt(n()))
    )

# Average expression heatmap
top_genes <- DEGs %>%
    filter(p_val_adj < 0.05, abs(avg_log2FC) > 1) %>%
    group_by(cluster) %>%
    slice_max(order_by = abs(avg_log2FC), n = 30) %>%
    pull(gene) %>%
    unique()

avg_exp <- AverageExpression(Metacells_2, features = top_genes,
                             group.by = "TwoClusters", assay = "RNA",
                             slot = "data")$RNA

col_fun <- colorRamp2(c(-2, 0, 2), c("navy", "white", "firebrick3"))
ht <- Heatmap(
    scale(t(avg_exp)),
    name            = "Scaled Avg Exp",
    col             = col_fun,
    cluster_rows    = TRUE,
    cluster_columns = TRUE,
    show_row_names  = TRUE,
    show_column_names = TRUE
)

pdf(file.path(output_dir, "Figure2_D.pdf"), width = 10, height = 6)
draw(ht)
dev.off()

# Write tables
DEGs$log2FC <- DEGs$avg_log2FC
write.table(DEGs, file.path(output_dir,
    "3.Table_DEGs_MTC_456_vs_123_FULL_nonsig_incl_tstat.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE)

genes_irea <- DEGs %>%
    filter(p_val_adj < 0.05, avg_log2FC > 1, cluster == "MTC_456") %>%
    select(gene, avg_log2FC, p_val_adj)
write.csv(genes_irea, file.path(output_dir, "2024_04_01_Genes_for_IREA.csv"),
          row.names = FALSE)

message("Saved outputs to: ", output_dir)
