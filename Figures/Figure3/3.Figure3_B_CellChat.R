#!/usr/bin/env Rscript
# Author: Dimitrios Kyriakis
# Description: CellChat cell-cell communication analysis on immune metacells.
#              Computes communication probabilities and plots circle diagrams
#              for COMPLEMENT, SPP1, and ApoE signalling pathways.
#              Produces Figure 3B circle plots.
# Usage: Rscript 3.Figure3_B_CellChat.R <metacell_rds> <output_dir>
#   <metacell_rds>  Path to 2024_03_28_Immune_Combined_Metacells.rds.ztsd
#   <output_dir>    Directory for output PDFs (default: current dir)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
    stop("Usage: Rscript 3.Figure3_B_CellChat.R <metacell_rds> [output_dir]")
}

input_rds  <- args[1]
output_dir <- if (length(args) >= 2) args[2] else "."

set.seed(12345)

suppressPackageStartupMessages({
    library(Seurat)
    library(CellChat)
    library(ggplot2)
    library(cowplot)
    library(ggplotify)
    library(patchwork)
})

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

Immune.combined <- readCRDS(input_rds)
Idents(Immune.combined) <- "Clusters"
Immune.combined$samples <- Immune.combined$orig.ident

# Node colours matching original analysis
node_colors <- c(
    "B Cells"   = "#F4B462",
    "MONO"      = "#80B1D3",
    "MTC_123"   = "#999999",
    "MTC_456"   = "#A32A31",
    "MTC_Prolif"= "#EE8172",
    "NEUT"      = "#DBDB98",
    "T Cells"   = "#B1DD69"
)

# Build CellChat object
cellchat <- createCellChat(object = Immune.combined,
                           group.by = "Clusters",
                           assay    = "RNAv3")

CellChatDB      <- CellChatDB.human
cellchat@DB     <- CellChatDB
cellchat        <- subsetData(cellchat)
cellchat        <- identifyOverExpressedGenes(cellchat)
cellchat        <- identifyOverExpressedInteractions(cellchat)
cellchat        <- computeCommunProb(cellchat, type = "triMean")
cellchat        <- filterCommunication(cellchat, min.cells = 10)
cellchat        <- computeCommunProbPathway(cellchat)
cellchat        <- aggregateNet(cellchat)
cellchat        <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")

groupSize <- as.numeric(table(cellchat@idents))

# Circle plots for three pathways
pathways_to_plot <- c("COMPLEMENT", "SPP1", "ApoE")
for (pw in pathways_to_plot) {
    out_file <- file.path(output_dir,
                          paste0("Figure3_C_", gsub("-", "_", pw), "_Pathway.pdf"))
    pdf(out_file, width = 8, height = 8)
    netVisual_aggregate(cellchat,
                        signaling          = pw,
                        layout             = "circle",
                        color.use          = node_colors,
                        vertex.label.cex   = 0)
    dev.off()
    message("Saved: ", out_file)
}

# Aggregate interaction weight diagram
pdf(file.path(output_dir, "Figure3_B_AllInteractions.pdf"), width = 7, height = 7)
groupSize <- as.numeric(table(cellchat@idents))
netVisual_circle(cellchat@net$weight,
                 vertex.weight = groupSize,
                 color.use     = node_colors,
                 weight.scale  = TRUE,
                 label.edge    = FALSE,
                 title.name    = "Interaction weights/strength")
dev.off()

message("Done. Outputs in: ", output_dir)
