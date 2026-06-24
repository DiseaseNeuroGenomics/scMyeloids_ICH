#!/usr/bin/env Rscript
# Author: Dimitrios Kyriakis
# Description: Compositional analysis of immune cell proportions using crumblr.
#              Tests mRS (good vs poor outcome) effect on cell-type composition.
#              Produces Figure 1D (tree + coefficient plots).
# Usage: Rscript 1.Figure1_D_CrumblR.R <seurat_rds> <threads> <output_pdf1> <output_pdf2>
#   <seurat_rds>   Path to 2024_11_All_Immune.rds
#   <threads>      Number of parallel threads (default: 4)
#   <output_pdf1>  Tree plot PDF (default: Figure1_E1.pdf)
#   <output_pdf2>  Coefficient plot PDF (default: Figure1_E2.pdf)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
    stop("Usage: Rscript 1.Figure1_D_CrumblR.R <seurat_rds> [threads] [output_pdf1] [output_pdf2]")
}

input_rds   <- args[1]
threads     <- if (length(args) >= 2) as.integer(args[2]) else 4L
output_pdf1 <- if (length(args) >= 3) args[3] else "Figure1_E1.pdf"
output_pdf2 <- if (length(args) >= 4) args[4] else "Figure1_E2.pdf"

set.seed(12345)

suppressPackageStartupMessages({
    library(Seurat)
    library(SingleCellExperiment)
    library(dreamlet)
    library(crumblr)
    library(variancePartition)
    library(muscat)
    library(ggtree)
    library(aplot)
    library(ggplot2)
    library(cowplot)
    library(BiocParallel)
    library(doParallel)
})

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

cobj <- crumblr(cellCounts(pbObj))

form <- ~ (1 | race) + (1 | sex) + (1 | BL) + (1 | dx) + age + TSH + mRS_binned + 0
L <- makeContrastsDream(
    form,
    colData(pbObj),
    contrasts = c(mRS_contrast = "mRS_binnedB - mRS_binnedG")
)

fit <- dream(cobj, form, colData(pbObj), L = L)
fit <- eBayes(fit)

res.mRS <- topTable(fit, coef = "mRS_contrast", number = Inf, sort.by = "none")
res.mRS$assay <- rownames(res.mRS)

res.meta <- dh_meta_analysis(list(res.mRS))
hc       <- buildClusterTreeFromPB(pbObj)
fig.tree <- dh_plotTree(ape::as.phylo(hc), xmax.scale = 2.2) +
    theme(legend.position = "bottom")
fig.es <- dh_plotCoef(res.meta, coef = "mRS", fig.tree, ylab = "mRS (Poor vs Good)")

pdf(output_pdf1, width = 8, height = 4)
print(fig.tree)
dev.off()

pdf(output_pdf2, width = 8, height = 4)
print(fig.es %>% insert_left(fig.tree, width = 1.4))
dev.off()

message("Saved: ", output_pdf1, " and ", output_pdf2)
