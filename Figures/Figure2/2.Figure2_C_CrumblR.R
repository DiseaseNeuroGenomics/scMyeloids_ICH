#!/usr/bin/env Rscript
# Author: Dimitrios Kyriakis
# Description: Compositional analysis of myeloid metacell subtypes using crumblr.
#              Tests mRS (good vs poor outcome) on myeloid subtype composition.
#              Produces Figure 2C (tree + coefficient plots).
# Usage: Rscript 2.Figure2_C_CrumblR.R <metacell_rds> <threads> <output_pdf1> <output_pdf2>
#   <metacell_rds>  Path to 2024_03_28_Myeloid_Metacells_Subclass_ADAM.rds.ztsd
#   <threads>       Number of parallel threads (default: 4)
#   <output_pdf1>   Tree plot PDF (default: Figure2_C1.pdf)
#   <output_pdf2>   Coefficient plot PDF (default: Figure2_C2.pdf)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
    stop("Usage: Rscript 2.Figure2_C_CrumblR.R <metacell_rds> [threads] [out1] [out2]")
}

input_rds   <- args[1]
threads     <- if (length(args) >= 2) as.integer(args[2]) else 4L
output_pdf1 <- if (length(args) >= 3) args[3] else "Figure2_C1.pdf"
output_pdf2 <- if (length(args) >= 4) args[4] else "Figure2_C2.pdf"

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

Metacells <- readCRDS(input_rds)

mat <- as.matrix(Metacells@assays$RNA@counts)
sce <- SingleCellExperiment(
    assays  = list(counts = mat),
    colData = Metacells@meta.data
)
rm(Metacells, mat)
gc()

pbObj <- aggregateToPseudoBulk(
    sce,
    assay      = "counts",
    cluster_id = "Subtype_Names",
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
hc        <- buildClusterTreeFromPB(pbObj)
fig.tree  <- dh_plotTree(ape::as.phylo(hc), xmax.scale = 2.2) +
    theme(legend.position = "bottom")
fig.es    <- dh_plotCoef(res.meta, coef = "mRS", fig.tree, ylab = "mRS (Poor vs Good)")

pdf(output_pdf1, width = 8, height = 4)
print(fig.tree)
dev.off()

pdf(output_pdf2, width = 8, height = 4)
print(fig.es %>% insert_left(fig.tree, width = 1.4))
dev.off()

message("Saved: ", output_pdf1, " and ", output_pdf2)
