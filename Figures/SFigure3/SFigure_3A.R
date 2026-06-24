#!/usr/bin/env Rscript
# Author: Dimitrios Kyriakis
# Description: Stacked bar plots showing myeloid subclass and subtype
#              composition per metacell cluster (dittoSeq).
#              Produces SFigure 3A.
# Usage: Rscript SFigure_3A.R <metacell_rds> <output_pdf>
#   <metacell_rds>  Path to 7.WGCNA/Metacells_Single_Level.rds.ztsd
#   <output_pdf>    Output PDF path (default: SFigure3A.pdf)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
    stop("Usage: Rscript SFigure_3A.R <metacell_rds> [output_pdf]")
}

input_rds  <- args[1]
output_pdf <- if (length(args) >= 2) args[2] else "SFigure3A.pdf"

set.seed(123456789)

suppressPackageStartupMessages({
    library(Seurat)
    library(dittoSeq)
    library(ggplot2)
    library(cowplot)
    library(patchwork)
    library(RColorBrewer)
})

source('utils.R')

Metacells_Single_Level <- readCRDS(input_rds)

p1 <- dittoBarPlot(
    Metacells_Single_Level,
    "Subclass_withoutMG",
    group.by    = "Metacells_identity",
    color.panel = Project_Colors
) +
    theme(text = element_text(size = 20)) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    ggtitle("Subclass composition per metacell cluster")

p2 <- dittoBarPlot(
    Metacells_Single_Level,
    "Subtype_withoutMG",
    group.by    = "Metacells_identity",
    color.panel = Project_Colors
) +
    theme(text = element_text(size = 20)) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    ggtitle("Subtype composition per metacell cluster")

pdf(output_pdf, width = 12, height = 5)
print(p1 + p2)
dev.off()

message("Saved: ", output_pdf)
