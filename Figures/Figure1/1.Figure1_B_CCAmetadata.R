#!/usr/bin/env Rscript
# Author: Dimitrios Kyriakis
# Description: Canonical correlation analysis of metadata variables across
#              cell populations. Produces Figure 1B heatmap (pheatmap).
# Usage: Rscript 1.Figure1_B_CCAmetadata.R <metadata_tsv> <output_pdf>
#   <metadata_tsv>  Path to 2024_11_metadata_new.tsv
#   <output_pdf>    Output PDF path (default: Figure1B_CCA.pdf)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
    stop("Usage: Rscript 1.Figure1_B_CCAmetadata.R <metadata_tsv> [output_pdf]")
}

input_tsv  <- args[1]
output_pdf <- if (length(args) >= 2) args[2] else "Figure1B_CCA.pdf"

set.seed(12345)

suppressPackageStartupMessages({
    library(variancePartition)
    library(pheatmap)
    library(RColorBrewer)
    library(dplyr)
    library(readr)
})

metadata <- as.data.frame(read_tsv(input_tsv))
rownames(metadata) <- metadata[[1]]
metadata[[1]] <- NULL

# Keep numeric / factor columns suitable for CCA
meta_num <- metadata[, sapply(metadata, function(x) is.numeric(x) | is.factor(x)), drop = FALSE]

# Canonical correlation between pairs of variables
canCor <- canCorPairs(~ ., meta_num)

# Plot
colors <- colorRampPalette(rev(brewer.pal(9, "RdBu")))(100)

pdf(output_pdf, width = 8, height = 7)
pheatmap(canCor,
         color        = colors,
         breaks       = seq(0, 1, length.out = 101),
         display_numbers = TRUE,
         number_format   = "%.2f",
         fontsize_number = 8,
         main         = "Canonical Correlation — Metadata Variables")
dev.off()

message("Saved: ", output_pdf)
