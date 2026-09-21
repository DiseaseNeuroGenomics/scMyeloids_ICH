# =============================================================================
# Figure 3C — CellChat pathway networks (Complement / SPP1 / ApoE)
# Reads:  a metacell Seurat object (.rds, or .ztsd/.zst compressed)
# Writes: Table17_CCI_CellChat.csv
#         Figures/Figure3_C_Complement_Pathway.pdf
#         Figures/Figure3_C_Spp1_Pathway.pdf
#         Figures/Figure3_C_Apoe_Pathway.pdf
# Run from inside the Liana/ folder.
# Requires the command-line `zstd` tool if the input is compressed:
#   brew install zstd
# =============================================================================
set.seed(12345)

library(Seurat)
library(CellChat)
library(dplyr)

# ----------------------------------------------------- 0. robust file reader ----
# Bypasses the R `archive` package (unreliable zstd support on some macOS
# libarchive builds) and shells out to the standalone zstd CLI instead.
read_metacells <- function(filename) {
    con <- archive::file_read(file = filename)
    res <- readRDS(con)
    close(con)
    res
}

RDS_IN  <- "2024_03_28_Immune_Combined_Metacells.rds.ztsd"  # set to whatever file is in this folder
OUT_DIR <- "Figures"
dir.create(OUT_DIR, showWarnings = FALSE)

node_colors <- c(
  "B Cells"    = "#F4B462",
  "MONO"       = "#80B1D3",
  "MTC_123"    = "#999999",
  "MTC_456"    = "#A32A31",
  "MTC_Prolif" = "#EE8172",
  "NEUT"       = "#DBDB98",
  "T Cells"    = "#B1DD69"
)

# ----------------------------------------------------- 1. build CellChat obj ----
Immune.combined <- read_metacells(RDS_IN)
Idents(Immune.combined) <- "Clusters"
Immune.combined$samples <- Immune.combined$orig.ident

cellchat <- createCellChat(object = Immune.combined, group.by = "Clusters")
cellchat@DB <- CellChatDB.human

cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

# ----------------------------------------------------- 2. infer communication ----
cellchat <- computeCommunProb(cellchat, type = "triMean")
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)

node_colors <- node_colors[colnames(cellchat@net$weight)]

# ----------------------------------------------------- 3. save the full table ----
TABLE_OUT <- "Table17_CCI_CellChat.csv"

full_table <- subsetCommunication(cellchat, thresh = 1)  # ALL source-target-LR triples, not filtered to significant
write.csv(full_table, TABLE_OUT, row.names = FALSE)
cat("Saved", nrow(full_table), "rows to", TABLE_OUT, "\n")

# ----------------------------------------------------- 4. save the 3 panels ----
pdf(file.path(OUT_DIR, "Figure3_C_Complement_Pathway.pdf"))
netVisual_aggregate(cellchat, signaling = "COMPLEMENT", layout = "circle",
                    color.use = node_colors, vertex.label.cex = 0,
                    edge.width.max = 5, arrow.size = 0.5, arrow.width = 1.8)
dev.off()

pdf(file.path(OUT_DIR, "Figure3_C_Spp1_Pathway.pdf"))
netVisual_aggregate(cellchat, signaling = "SPP1", layout = "circle",
                    color.use = node_colors, vertex.label.cex = 0,
                    edge.width.max = 5, arrow.size = 0.5, arrow.width = 1.8)
dev.off()

pdf(file.path(OUT_DIR, "Figure3_C_Apoe_Pathway.pdf"))
netVisual_aggregate(cellchat, signaling = "ApoE", layout = "circle",
                    color.use = node_colors, vertex.label.cex = 0,
                    edge.width.max = 5, arrow.size = 0.5, arrow.width = 1.8)
dev.off()

cat("Saved 3 panels to", OUT_DIR, "\n")


# ----------------------------------------------------- 5. combined figure -----
# Complement / Spp1 / ApoE side-by-side, one shared cell-type legend.
COMBINED_OUT <- file.path(OUT_DIR, "Figure3_C_Combined.pdf")

pathways_to_plot <- c("COMPLEMENT", "SPP1", "ApoE")
panel_titles     <- c("Complement", "Spp1", "ApoE")

pdf(COMBINED_OUT, width = 13, height = 4.5)
layout(matrix(1:4, nrow = 1), widths = c(1, 1, 1, 0.45))
par(mar = c(1, 1, 3, 1), xpd = TRUE)

for (i in seq_along(pathways_to_plot)) {
    netVisual_aggregate(cellchat, signaling = pathways_to_plot[i], layout = "circle",
                        color.use = node_colors, vertex.label.cex = 0,
                        edge.width.max = 5, arrow.size = 0.5, arrow.width = 1.8)
    title(panel_titles[i], line = 0.5)
  }

plot.new()
legend("left", legend = names(node_colors), pch = 21,
      pt.bg = node_colors, pt.cex = 2, bty = "n", cex = 0.9,
      title = "Cell type")

dev.off()
cat("Saved combined panel to", COMBINED_OUT, "\n")
