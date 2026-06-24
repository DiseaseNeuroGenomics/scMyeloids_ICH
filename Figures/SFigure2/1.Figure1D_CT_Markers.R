#!/usr/bin/env Rscript
# Author: Dimitrios Kyriakis
# Description: ComplexHeatmap of top cell-type marker genes across major immune
#              cell populations. Uses cellTypeSpecificity to rank genes.
#              Produces Figure 1D.
# Usage: Rscript 1.Figure1D_CT_Markers.R <seurat_rds> <output_pdf>
#   <seurat_rds>   Path to 2024_11_All_Immune.rds
#   <output_pdf>   Output PDF path (default: Figure1_D.pdf)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
    stop("Usage: Rscript 1.Figure1D_CT_Markers.R <seurat_rds> [output_pdf]")
}

input_rds  <- args[1]
output_pdf <- if (length(args) >= 2) args[2] else "Figure1_D.pdf"

set.seed(12345)

suppressPackageStartupMessages({
    library(Seurat)
    library(SingleCellExperiment)
    library(dreamlet)
    library(muscat)
    library(ComplexHeatmap)
    library(circlize)
    library(RColorBrewer)
    library(BiocParallel)
})

Immune <- readRDS(input_rds)

mat <- as.matrix(Immune@assays$RNA@counts)
sce <- SingleCellExperiment(
    assays  = list(counts = mat),
    colData = Immune@meta.data
)
rm(mat); gc()

pbObj_Immune <- aggregateToPseudoBulk(
    sce,
    assay      = "counts",
    cluster_id = "CellType",
    sample_id  = "donor",
    BPPARAM    = SnowParam(4, progressbar = TRUE)
)

n_genes <- 10
df_cts  <- cellTypeSpecificity(pbObj_Immune)
df_cts  <- df_cts[df_cts$totalCPM > 100, ]
df_cts  <- df_cts[grep("\\.", rownames(df_cts), invert = TRUE), ]

top_idx <- function(x, n) order(x, decreasing = TRUE)[seq_len(n)]
genes   <- rownames(df_cts)[apply(df_cts, 2, top_idx, n = n_genes)]
genes   <- genes[-seq_len(n_genes)]

avg_ex <- AverageExpression(Immune, group.by = "CellType", return.seurat = TRUE)
mat_sc <- t(avg_ex@assays$RNA@scale.data[genes, , drop = FALSE])

split <- rep(rownames(mat_sc), each = n_genes)

bar_colors <- c("#FDB462", "#80B1D3", "#FB8072", "#B3DE69", "#8DD3C7")
labels_h   <- c("B", "MONO", "Myeloid", "NEUT", "TCells")

ha <- HeatmapAnnotation(
    foo = anno_block(gp     = gpar(fill = bar_colors),
                     labels = labels_h)
)

ha1 <- rowAnnotation(
    CL = factor(rownames(mat_sc)),
    col = list(CL = setNames(bar_colors, labels_h))
)

col_fun <- colorRamp2(c(-3, 0, 3), c("blue", "white", "red"))

ht <- Heatmap(
    mat_sc,
    name              = "Avg. exp.",
    col               = col_fun,
    cluster_columns   = FALSE,
    cluster_rows      = FALSE,
    column_split      = split,
    column_order      = genes,
    top_annotation    = ha,
    right_annotation  = ha1,
    row_title         = NULL,
    column_title      = NULL,
    heatmap_legend_param = list(
        at              = c(-4, -2, 0, 2, 4),
        title           = "Avg. exp.",
        legend_height   = unit(4, "cm"),
        title_position  = "leftcenter-rot"
    )
)

pdf(output_pdf, width = 12, height = 6)
draw(ht)
dev.off()

message("Saved: ", output_pdf)
