#!/usr/bin/env Rscript
# Author: Dimitrios Kyriakis
# Description: SCENIC regulon activity analysis. Loads pySCENIC loom output,
#              extracts TF AUC scores, runs Fisher exact tests between regulon
#              targets and DEG lists, and plots TF log2FC vs odds-ratio scatter.
#              Produces Figure 2E.
# Usage: Rscript 2.Figure2_E_GRN.R <loom_file> <tf_csv> <deg_csv> <output_pdf>
#   <loom_file>   Path to SCENIC output loom file (SCopeLoomR format)
#   <tf_csv>      Path to 2025_SCENIC_DEG_TFs.csv
#   <deg_csv>     Path to 2025_SCENIC_DEGs_genes.csv
#   <output_pdf>  Output PDF path (default: Figure2_F.pdf)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
    stop("Usage: Rscript 2.Figure2_E_GRN.R <loom_file> <tf_csv> <deg_csv> [output_pdf]")
}

loom_file  <- args[1]
tf_csv     <- args[2]
deg_csv    <- args[3]
output_pdf <- if (length(args) >= 4) args[4] else "Figure2_F.pdf"

set.seed(12345)

suppressPackageStartupMessages({
    library(SCopeLoomR)
    library(ggplot2)
    library(cowplot)
    library(dplyr)
    library(RColorBrewer)
})

loom <- open_loom(loom_file)

regulons   <- get_regulons(loom, column.attr.name = "Regulons")
auc_mtx    <- get_regulons_AUC(loom, column.attr.name = "RegulonsAUC")
regulon_df <- as.data.frame(regulons)
close_loom(loom)

# Load pre-computed TF differential stats and DEG background
tf_df  <- read.csv(tf_csv, stringsAsFactors = FALSE)
deg_df <- read.csv(deg_csv, stringsAsFactors = FALSE)

background <- unique(deg_df$gene)

# Fisher exact test: for each TF regulon, test enrichment in DEG lists
fisher_results <- lapply(rownames(regulon_df), function(tf) {
    targets <- names(which(regulon_df[tf, ] == 1))
    for_each <- lapply(c("MTC_456", "MTC_123"), function(cl) {
        deg_cl <- deg_df$gene[deg_df$cluster == cl & deg_df$p_val_adj < 0.05]
        a <- length(intersect(targets, deg_cl))
        b <- length(setdiff(targets, deg_cl))
        cc <- length(setdiff(deg_cl, targets))
        d <- length(setdiff(background, union(targets, deg_cl)))
        ft <- fisher.test(matrix(c(a, b, cc, d), nrow = 2))
        data.frame(TF = tf, cluster = cl,
                   OddsRatio = ft$estimate,
                   p.value   = ft$p.value,
                   stringsAsFactors = FALSE)
    })
    do.call(rbind, for_each)
})
fisher_df <- do.call(rbind, fisher_results)
fisher_df$FDR <- p.adjust(fisher_df$p.value, method = "BH")

# Merge with TF log2FC
merged <- merge(fisher_df, tf_df[, c("gene", "avg_log2FC", "cluster")],
                by.x = c("TF", "cluster"), by.y = c("gene", "cluster"))
merged$Significant <- merged$FDR < 0.05

p <- ggplot(merged, aes(x = avg_log2FC, y = log2(OddsRatio + 1e-6),
                         colour = cluster, size = -log10(FDR + 1e-300))) +
    geom_point(alpha = 0.7) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_vline(xintercept = 0, linetype = "dashed") +
    scale_colour_manual(values = c("MTC_123" = "#4575B4", "MTC_456" = "#D73027")) +
    labs(x = "TF avg log2FC (MTC_456 vs MTC_123)",
         y = "log2(Odds Ratio) — regulon target enrichment",
         colour = "MTC Group",
         size   = "-log10(FDR)") +
    theme_cowplot() +
    theme(legend.position = "right")

pdf(output_pdf, width = 8, height = 6)
print(p)
dev.off()

message("Saved: ", output_pdf)
