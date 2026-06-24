#!/usr/bin/env Rscript
# Author: Dimitrios Kyriakis
# Description: UMAP DimPlot of all immune cells coloured by cell type.
#              Produces Figure 1C.
# Usage: Rscript 1.Figure1_C_UMAP.R <seurat_rds> <output_pdf>
#   <seurat_rds>   Path to 2024_11_All_Immune.rds
#   <output_pdf>   Output PDF path (default: Figure1C_UMAP.pdf)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
    stop("Usage: Rscript 1.Figure1_C_UMAP.R <seurat_rds> [output_pdf]")
}

input_rds  <- args[1]
output_pdf <- if (length(args) >= 2) args[2] else "Figure1C_UMAP.pdf"

set.seed(12345)

suppressPackageStartupMessages({
    library(Seurat)
    library(scCustomize)
    library(RColorBrewer)
    library(ggplot2)
    library(cowplot)
})

source('utils.R')

Immune <- readRDS(input_rds)

p <- DimPlot_scCustom(
    seurat_object = Immune,
    pt.size       = 0.05,
    colors_use    = Project_Colors,
    group.by      = "CellType",
    reduction     = "umap",
    figure_plot   = TRUE
)

pdf(output_pdf, width = 10, height = 8)
print(p)
dev.off()

message("Saved: ", output_pdf)
