#!/usr/bin/env Rscript
# Author: Dimitrios Kyriakis
# Description: Three-panel SFigure 2:
#   B — Marker gene heatmap for FreshMG-mapped myeloid subtypes.
#   C — crumblr compositional analysis (mRS) of myeloid subtypes.
#   D — Dreamlet DE bubble plot across myeloid subtypes.
# Usage: Rscript SFigure2_BCD.R <seurat_rds> <threads> <output_dir>
#   <seurat_rds>   Path to 2024_11_All_Immune.rds
#   <threads>      Number of parallel threads (default: 4)
#   <output_dir>   Output directory (default: SupFigures/)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
    stop("Usage: Rscript SFigure2_BCD.R <seurat_rds> [threads] [output_dir]")
}

input_rds  <- args[1]
threads    <- if (length(args) >= 2) as.integer(args[2]) else 4L
output_dir <- if (length(args) >= 3) args[3] else "SupFigures"

set.seed(12345)

suppressPackageStartupMessages({
    library(Seurat)
    library(SingleCellExperiment)
    library(dreamlet)
    library(crumblr)
    library(variancePartition)
    library(muscat)
    library(ComplexHeatmap)
    library(circlize)
    library(ggtree)
    library(aplot)
    library(ggplot2)
    library(cowplot)
    library(dplyr)
    library(RColorBrewer)
    library(BiocParallel)
    library(doParallel)
})

source('utils.R')

registerDoParallel(cores = threads)
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

Immune <- readRDS(input_rds)

# Subset myeloid only
Microglia_only <- subset(Immune, subset = CellPop %in% c("Myeloid"))
rm(Immune); gc()

# ── Panel B: Marker gene heatmap ──────────────────────────────────────────────
my_vector <- c(
    "AIF1", "CCL3", "CCL4", "HIF1A", "HSPA1A", "HSPH1", "HSP90AA1", "DNAJA4",
    "IFIT1", "IFIT2", "IFIT3", "TMEM163", "NCK2", "SH3RF3",
    "FRMD4A", "P2RY12", "CX3CR1", "PICALM",
    "CD163", "CSGALNACT1", "GPNMB", "ERN1", "CSKMT", "SAMD4A",
    "MKI67", "HELLS", "CLSPN", "CENPK", "PLK2"
)

ToPlot <- AverageExpression(
    Microglia_only, group.by = "Subtype", assay = "RNA",
    slot = "data", features = intersect(my_vector, rownames(Microglia_only))
)$RNA

subtype_order <- c(
    "Adapt_AIF1", "Adapt_CCL3", "Adapt_HIF1A", "Adapt_HIST",
    "Adapt_HSPA1A", "Adapt_IFI44L", "Adapt_TMEM163",
    "Homeo_FRMD4A", "Homeo_PICALM", "PVM_CD163", "ADAM_GPNMB",
    "Prolif_MKI67", "exAM_ERN1"
)
subtype_order <- intersect(subtype_order, colnames(ToPlot))

z_scores <- scale(t(ToPlot[, subtype_order, drop = FALSE]))

col_fun <- colorRamp2(c(-2, 0, 2), c("navy", "white", "firebrick3"))
ht <- Heatmap(z_scores,
    name              = "Mean Z-Score",
    col               = col_fun,
    cluster_rows      = FALSE,
    cluster_columns   = FALSE,
    show_row_names    = TRUE,
    show_column_names = TRUE,
    row_names_side    = "left",
    column_names_rot  = 90
)

pdf(file.path(output_dir, "SFigure2B.pdf"), width = 15, height = 12)
draw(ht)
dev.off()
message("Saved SFigure2B")

# ── Panel C: crumblr compositional analysis ───────────────────────────────────
Microglia_only$Subtype_Names <- gsub(
    "Adapt_|Homeo_|PVM_|ADAM_|Prolif_|exAM_", "", Microglia_only$Subtype
)

mat <- as.matrix(Microglia_only@assays$RNA@counts)
sce <- SingleCellExperiment(
    assays  = list(counts = mat),
    colData = Microglia_only@meta.data
)
rm(mat); gc()

pbObj_subtype <- aggregateToPseudoBulk(
    sce, assay = "counts",
    cluster_id = "Subtype_Names", sample_id = "donor",
    BPPARAM    = SnowParam(threads, progressbar = TRUE)
)

cobj <- crumblr(cellCounts(pbObj_subtype))

form <- ~ (1 | race) + (1 | sex) + (1 | BL) + (1 | dx) + age + TSH + mRS_binned + 0
L <- makeContrastsDream(
    form, colData(pbObj_subtype),
    contrasts = c(mRS_contrast = "mRS_binnedB - mRS_binnedG")
)
fit  <- dream(cobj, form, colData(pbObj_subtype), L = L)
fit  <- eBayes(fit)
res.mRS <- topTable(fit, coef = "mRS_contrast", number = Inf, sort.by = "none")
res.mRS$assay <- rownames(res.mRS)

res.meta <- dh_meta_analysis(list(res.mRS))
hc       <- buildClusterTreeFromPB(pbObj_subtype)
fig.tree <- dh_plotTree(ape::as.phylo(hc), xmax.scale = 2.2) +
    theme(legend.position = "bottom")
fig.es   <- dh_plotCoef(res.meta, coef = "mRS", fig.tree, ylab = "mRS")

pdf(file.path(output_dir, "SFigure2C.pdf"), width = 8, height = 4)
print(fig.es %>% insert_left(fig.tree, width = 1.4))
dev.off()
message("Saved SFigure2C")

# ── Panel D: Dreamlet bubble plot ─────────────────────────────────────────────
pbObj_sub2 <- aggregateToPseudoBulk(
    sce, assay = "counts",
    cluster_id = "Subtype", sample_id = "donor",
    BPPARAM    = SnowParam(threads, progressbar = TRUE)
)

res_proc <- processAssays(
    pbObj_sub2,
    ~ race + dx + BL + TSH + age + mRS_binned,
    min.count = 1,
    BPPARAM   = SnowParam(threads, progressbar = TRUE)
)

res_dream <- dreamlet(
    res_proc,
    ~ (1 | race) + (1 | dx) + (1 | BL) + age + mRS_binned + 0,
    contrasts = c(Diff_B_vs_G = "mRS_binnedB-mRS_binnedG"),
    BPPARAM   = SnowParam(threads, progressbar = TRUE)
)

DEGs_file <- as.data.frame(
    topTable(res_dream, p.value = 1, lfc = 0, number = Inf, coef = "Diff_B_vs_G")
)
DEGs_file <- DEGs_file %>%
    group_by(ID) %>%
    mutate(avg_LFC    = sum(logFC) / length(unique(DEGs_file$assay)),
           pos_avgLFC = mean(pmax(logFC, 0)))

sign_df   <- DEGs_file[DEGs_file$adj.P.Val < 0.05 & abs(DEGs_file$avg_LFC) > 1, ]
data_plot <- DEGs_file[DEGs_file$ID %in% unique(sign_df$ID), ] %>% na.omit()
data_plot <- data_plot[order(data_plot$avg_LFC), ]
data_plot$gene  <- factor(data_plot$ID, levels = unique(data_plot$ID))
data_plot$assay <- factor(data_plot$assay,
    levels = rev(subtype_order[subtype_order %in% unique(data_plot$assay)])
)

bubble <- ggplot(data_plot,
    aes(x = gene, y = assay, size = abs(logFC), fill = logFC)) +
    geom_point(shape = 21, colour = "black") +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red") +
    theme_cowplot() +
    theme(
        panel.border    = element_rect(colour = "black", fill = NA, linewidth = 1),
        axis.text.x     = element_text(angle = 90, vjust = 0.5, hjust = 1),
        legend.position = "bottom"
    ) +
    labs(x = NULL, y = NULL)

pdf(file.path(output_dir, "SFigure2D.pdf"), width = 10, height = 5)
print(bubble)
dev.off()
message("Saved SFigure2D. All done.")
