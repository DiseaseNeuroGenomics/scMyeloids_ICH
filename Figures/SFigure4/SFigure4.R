#!/usr/bin/env Rscript
# Author: Dimitrios Kyriakis
# Description: Supplementary Figure 4 — label transfer from mouse MCAO scRNA-seq
#              (GSE189432) to ICH myeloid metacells via Seurat FindTransferAnchors.
# Usage: Rscript SFigure4.R <external_data_dir> <metacells_rds> <output_dir>
#   <external_data_dir>  Directory containing GSE189432 files and annotation CSV
#   <metacells_rds>      Path to myeloid metacell RDS (.rds.ztsd)
#   <output_dir>         Directory for output PDF/PNG

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
    stop("Usage: Rscript SFigure4.R <external_data_dir> <metacells_rds> <output_dir>")
}

external_data_dir <- args[1]
metacells_rds     <- args[2]
output_dir        <- args[3]

suppressPackageStartupMessages({
    library(GEOquery)
    library(Seurat)
    library(tibble)
    library(magrittr)
    library(dittoSeq)
    library(patchwork)
    library(harmony)
})

source('scripts/utils.R')
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

set.seed(24022012)

# ===================== Load GSE189432 annotation =====================
Metadata <- read.csv(file.path(external_data_dir, "GSE189432_annotations.csv.gz"))
Metadata <- Metadata[Metadata$sample %in% c('ctrl_cns','stroke_cns_24h_1','stroke_cns_24h_2','stroke_cns_72h'), ]

# ===================== Read 10X data =====================
read_sample <- function(name, mtx, feat, bar, suffix) {
    obj <- CreateSeuratObject(ReadMtx(
        mtx      = file.path(external_data_dir, mtx),
        features = file.path(external_data_dir, feat),
        cells    = file.path(external_data_dir, bar)
    ))
    obj$orig.ident <- name
    new_barcodes   <- gsub("-1", suffix, rownames(obj@meta.data))
    RenameCells(obj, new.names = new_barcodes)
}

CTRL  <- read_sample('ctrl_cns',
    'ctrl_cns/GSM5701739_ctrl_cns_matrix.mtx.gz',
    'ctrl_cns/GSM5701739_ctrl_cns_features.tsv.gz',
    'ctrl_cns/GSM5701739_ctrl_cns_barcodes.tsv.gz',
    '_ctrl_cns')
S24_1 <- read_sample('stroke_cns_24h_1',
    'stroke_24_cns/GSM5701742_stroke_cns_24h_1_matrix.mtx.gz',
    'stroke_24_cns/GSM5701742_stroke_cns_24h_1_features.tsv.gz',
    'stroke_24_cns/GSM5701742_stroke_cns_24h_1_barcodes.tsv.gz',
    '_stroke_cns_24h_1')
S24_2 <- read_sample('stroke_cns_24h_2',
    'stroke_24_cns_2/GSM5701743_stroke_cns_24h_2_matrix.mtx.gz',
    'stroke_24_cns_2/GSM5701743_stroke_cns_24h_2_features.tsv.gz',
    'stroke_24_cns_2/GSM5701743_stroke_cns_24h_2_barcodes.tsv.gz',
    '_stroke_cns_24h_2')
S72   <- read_sample('stroke_cns_72h',
    'stroke_72_cns/GSM5701746_stroke_cns_72h_matrix.mtx.gz',
    'stroke_72_cns/GSM5701746_stroke_cns_72h_features.tsv.gz',
    'stroke_72_cns/GSM5701746_stroke_cns_72h_barcodes.tsv.gz',
    '_stroke_cns_72h')

MICE <- merge(x = CTRL, y = c(S24_1, S24_2, S72), project = "Mice")
MICE$barcode <- rownames(MICE@meta.data)
MICE <- subset(MICE, subset = barcode %in% Metadata$barcode)
MICE$UMAP_1   <- Metadata$UMAP_1
MICE$UMAP_2   <- Metadata$UMAP_2
MICE$cluster  <- Metadata$cluster
MICE$sample   <- Metadata$sample
MICE[["RNA"]] <- JoinLayers(MICE[["RNA"]])

# Uppercase gene names to match human convention
counts <- MICE@assays$RNA$counts
rownames(counts) <- toupper(rownames(MICE))
MICE   <- CreateSeuratObject(counts = counts, meta.data = MICE@meta.data)

umap_coords <- Metadata[, c('UMAP_1','UMAP_2')]
rownames(umap_coords) <- Metadata$barcode
MICE[["ref.umap"]] <- CreateDimReducObject(
    embeddings = as.matrix(umap_coords),
    key = "UMAP_", assay = DefaultAssay(MICE)
)

# ===================== SCTransform + Harmony =====================
ifnb <- MICE
ifnb[["RNA"]] <- JoinLayers(ifnb[["RNA"]])
ifnb <- SCTransform(ifnb)
ifnb <- RunPCA(ifnb)
ifnb <- RunHarmony(ifnb, group.by.vars = "sample", dims.use = 1:40)
ifnb <- RunUMAP(ifnb, dims = 1:40, reduction = "harmony", min.dist = 0.3,
                seed.use = 24022012, reduction.name = "umap.harmony",
                reduction.key = "Uh_", return.model = TRUE)

colors_paper <- c(
    "Micro_1"="#FFC312","Micro_2"="#C4E538","Micro_3"="#12CBC4",
    "stress_Micro"="#FDA7DF","CAM_1"="#ED4C67","CAM_2"="#F79F1F",
    "SAMC"="#A3CB38","Macro_1"="#1289A7","Macro_2"="#D980FA",
    "stress_Myeloid"="#B53471","mDC1"="#EE5A24","mDC2"="#009432",
    "Granulo_1"="#0652DD","Granulo_2"="#9980FA","Mast"="#833471",
    "prolif_cells"="#EA2027","Bc"="#006266","gdTc"="#1B1464",
    "ILC2"="#5758BB","Tc"="#6F1E51","NK"="#40407A"
)

# ===================== Label transfer to metacells =====================
Metacells <- readCRDS(metacells_rds)
Idents(Metacells) <- 'Clusters'

anchors <- FindTransferAnchors(
    reference           = ifnb,
    query               = Metacells,
    normalization.method = "SCT",
    reference.reduction = "pca",
    dims                = 1:40
)

Metacells <- MapQuery(
    anchorset          = anchors,
    query              = Metacells,
    reference          = ifnb,
    refdata            = list(celltype = "cluster"),
    reference.reduction = "harmony",
    reduction.model    = "umap.harmony"
)

# ===================== Plots =====================
p1 <- DimPlot(ifnb, group.by = 'cluster', reduction = 'umap.harmony',
              pt.size = 0.0001, cols = colors_paper, label = FALSE) +
    xlim(-15, 10) + ylim(-15, 10) + ggtitle('Reference')
p2 <- DimPlot(Metacells, reduction = "umap", group.by = "predicted.celltype",
              cols = colors_paper) + NoLegend()
p3 <- dittoBarPlot(Metacells, group.by = 'Clusters', var = 'predicted.celltype',
                   color.panel = colors_paper) + coord_flip() + NoLegend()

pdf(file.path(output_dir, "SupFig.4.pdf"), width = 18, height = 6)
print(p1 + p2 + p3 + plot_layout(widths = c(6, 6, 3), guides = 'collect'))
dev.off()

png(file.path(output_dir, "SupFig.4.png"), width = 1800, height = 600, res = 100)
print(p1 + p2 + p3 + plot_layout(widths = c(6, 6, 3), guides = 'collect'))
dev.off()

message("Saved outputs to: ", output_dir)
