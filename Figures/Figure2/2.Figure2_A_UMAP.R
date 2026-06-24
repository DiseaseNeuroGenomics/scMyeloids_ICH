#!/usr/bin/env Rscript
# Author: Dimitrios Kyriakis
# Description: UMAP DimPlot of myeloid metacells coloured by Subclass and
#              two-cluster assignment (MTC_123 / MTC_456). Produces Figure 2A.
# Usage: Rscript 2.Figure2_A_UMAP.R <metacell_rds> <output_pdf>
#   <metacell_rds>  Path to 2024_03_28_Myeloid_Metacells_Subclass_ADAM.rds.ztsd
#   <output_pdf>    Output PDF path (default: Figure2_B.pdf)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
    stop("Usage: Rscript 2.Figure2_A_UMAP.R <metacell_rds> [output_pdf]")
}

input_rds  <- args[1]
output_pdf <- if (length(args) >= 2) args[2] else "Figure2_B.pdf"

set.seed(12345)

suppressPackageStartupMessages({
    library(Seurat)
    library(scCustomize)
    library(RColorBrewer)
    library(ggplot2)
    library(cowplot)
    library(patchwork)
})

source('utils.R')

Metacells <- readCRDS(input_rds)

Metacells$TwoClusters <- as.vector(Metacells$Clusters)
Metacells$TwoClusters[Metacells$Clusters %in% c("MTC_1", "MTC_2", "MTC_3")] <- "MTC_123"
Metacells$TwoClusters[Metacells$Clusters %in% c("MTC_4", "MTC_5", "MTC_6")] <- "MTC_456"

p1 <- DimPlot(Metacells, group.by = "Subclass", pt.size = 0.5) +
    ggtitle("Myeloid Metacells — Subclass")

p2 <- DimPlot(Metacells, group.by = "TwoClusters",
              cols = c("MTC_123" = "#4575B4", "MTC_456" = "#D73027"),
              pt.size = 0.5) +
    ggtitle("Myeloid Metacells — MTC Groups")

pdf(output_pdf, width = 14, height = 6)
print(p1 + p2)
dev.off()

message("Saved: ", output_pdf)
