#!/usr/bin/env Rscript
# Author: Dimitrios Kyriakis
# Description: Pseudobulk differential expression analysis with dreamlet across
#              major immune cell populations (ICH good vs poor outcome).
#              Produces Figure 1E (volcano plots + bubble plot).
# Usage: Rscript 1.Figure1_E_Dreamlet.R <seurat_rds> <threads> <output_pdf>
#   <seurat_rds>   Path to 2024_11_All_Immune.rds
#   <threads>      Number of parallel threads (default: 4)
#   <output_pdf>   Output PDF path (default: Figure1_F.pdf)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
    stop("Usage: Rscript 1.Figure1_E_Dreamlet.R <seurat_rds> [threads] [output_pdf]")
}

input_rds  <- args[1]
threads    <- if (length(args) >= 2) as.integer(args[2]) else 4L
output_pdf <- if (length(args) >= 3) args[3] else "Figure1_F.pdf"

set.seed(12345)

suppressPackageStartupMessages({
    library(Seurat)
    library(SingleCellExperiment)
    library(dreamlet)
    library(variancePartition)
    library(muscat)
    library(ggplot2)
    library(cowplot)
    library(ggrepel)
    library(dplyr)
    library(BiocParallel)
    library(doParallel)
})

source('utils.R')

registerDoParallel(cores = threads)

Immune <- readRDS(input_rds)

mat <- as.matrix(Immune@assays$RNA@counts)
sce <- SingleCellExperiment(
    assays  = list(counts = mat),
    colData = Immune@meta.data
)
rm(Immune, mat)
gc()

pbObj <- aggregateToPseudoBulk(
    sce,
    assay      = "counts",
    cluster_id = "CellType",
    sample_id  = "donor",
    BPPARAM    = SnowParam(threads, progressbar = TRUE)
)

form <- ~ race + dx + BL + TSH + age + mRS_binned
res.proc <- processAssays(pbObj, form, min.count = 1,
                          BPPARAM = SnowParam(threads, progressbar = TRUE))

vp.lst <- fitVarPart(res.proc, ~ race + dx + sex + BL + TSH + age + mRS_binned)

res.dream <- dreamlet(
    res.proc,
    ~ (1 | race) + (1 | dx) + (1 | BL) + age + mRS_binned + 0,
    contrasts  = c(Diff_B_vs_G = "mRS_binnedB-mRS_binnedG"),
    BPPARAM    = SnowParam(threads, progressbar = TRUE)
)

volcano_df <- as.data.frame(
    topTable(res.dream, p.value = 1, lfc = 0, number = Inf, coef = "Diff_B_vs_G")
)
volcano_df <- volcano_df[order(volcano_df$adj.P.Val), ]
volcano_df$diffexpressed <- "NS"
volcano_df$diffexpressed[volcano_df$logFC > 0.5 & volcano_df$adj.P.Val < 0.05]  <- "UP"
volcano_df$diffexpressed[volcano_df$logFC < -0.5 & volcano_df$adj.P.Val < 0.05] <- "DOWN"

# Volcano plots per cell type
x_lim <- max(abs(volcano_df$logFC), na.rm = TRUE) * 1.05
plotlist <- lapply(unique(volcano_df$assay), function(ct) {
    my_Vlc_plot(volcano_df, x_lim, ct)
})

# Bubble plot of significant genes
volcano_df <- volcano_df %>%
    group_by(ID) %>%
    mutate(pos_avgLFC = mean(logFC[logFC > 0], na.rm = TRUE)) %>%
    mutate(avg_LFC    = sum(logFC) / length(unique(volcano_df$assay)))

sign_df   <- volcano_df[volcano_df$adj.P.Val < 0.05 & abs(volcano_df$avg_LFC) > 1, ]
data_plot <- volcano_df[volcano_df$ID %in% unique(sign_df$ID), ] %>% na.omit()
data_plot <- data_plot[order(data_plot$avg_LFC), ]
data_plot$gene <- factor(data_plot$ID, levels = unique(data_plot$ID))

bubble <- ggplot(data_plot, aes(x = gene, y = assay, size = abs(logFC), fill = logFC)) +
    geom_point(shape = 21, colour = "black") +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red") +
    theme_cowplot() +
    theme(
        panel.border    = element_rect(colour = "black", fill = NA, linewidth = 1),
        axis.text.x     = element_text(angle = 90, vjust = 0.5, hjust = 1),
        legend.position = "bottom"
    ) +
    labs(x = NULL, y = NULL)

pdf(output_pdf, width = 20, height = 15)
gridExtra::grid.arrange(grobs = plotlist, ncol = 5)
print(bubble)
dev.off()

message("Saved: ", output_pdf)
