# ============================================================================
# Figure2_F_00_build_SAMC_reference.R
# Project: Decoding Immune Landscapes and Therapeutic Targets in Intracerebral Hemorrhage
# Author:  Dimitrios Kyriakis
# ============================================================================

# ==============================================================================
# Figure 2, Panel F (input build) — construct the harmonized mouse stroke
# reference (Beuker et al. 2022, GSE189432), which panel F maps the human MTC
# metacells onto to test whether mouse stroke-associated myeloid cells (SAMC)
# align preferentially with MTC_456.
#
# NOTE on paths — anchored under BASE_DIR/External_Data. Confirmed flat layout
# (no per-sample subfolders): GSM-prefixed *_barcodes/features/matrix.gz files
# sit directly in External_Data/, one triplet per sample.
#
# Input:  BASE_DIR/External_Data/GSE189432_annotations.csv.gz
#         BASE_DIR/External_Data/GSM*_{barcodes,features,matrix}.{tsv,mtx}.gz
# Output: DERIVED_DIR/Figure2_F_SAMC_reference.rds
# ==============================================================================

source("../../00_setup.R")
library(harmony)

# JoinLayers() only applies to Seurat v5 "Assay5" objects. This environment's
# default (per Seurat's own startup message) is the legacy "Assay" object,
# where merge() already concatenates directly and there's nothing to join —
# calling JoinLayers() there errors with "no applicable method". Guard it so
# the script works under either object version.
join_layers_if_needed <- function(obj, assay = "RNA") {
  if (inherits(obj[[assay]], "Assay5")) obj[[assay]] <- JoinLayers(obj[[assay]])
  obj
}

EXTERNAL_DATA_DIR <- file.path(BASE_DIR, "External_Data")

Metadata <- read.csv(file.path(EXTERNAL_DATA_DIR, "GSE189432_annotations.csv.gz"))
Metadata <- Metadata[Metadata$sample %in% c("ctrl_cns", "stroke_cns_24h_1", "stroke_cns_24h_2", "stroke_cns_72h"), ]

read_sample <- function(gsm_prefix, orig_ident) {
  obj <- CreateSeuratObject(ReadMtx(
    mtx      = file.path(EXTERNAL_DATA_DIR, paste0(gsm_prefix, "_matrix.mtx.gz")),
    features = file.path(EXTERNAL_DATA_DIR, paste0(gsm_prefix, "_features.tsv.gz")),
    cells    = file.path(EXTERNAL_DATA_DIR, paste0(gsm_prefix, "_barcodes.tsv.gz"))
  ))
  obj$orig.ident <- orig_ident
  new_names <- gsub("-1", paste0("_", orig_ident), rownames(obj@meta.data))
  RenameCells(obj, new.names = new_names)
}

CTRL  <- read_sample("GSM5701739_ctrl_cns",         "ctrl_cns")
S24_1 <- read_sample("GSM5701742_stroke_cns_24h_1", "stroke_cns_24h_1")
S24_2 <- read_sample("GSM5701743_stroke_cns_24h_2", "stroke_cns_24h_2")
S72   <- read_sample("GSM5701746_stroke_cns_72h",   "stroke_cns_72h")

MICE <- merge(x = CTRL, y = list(S24_1, S24_2, S72), project = "Mice")
MICE$barcode <- rownames(MICE@meta.data)
MICE <- subset(MICE, subset = barcode %in% Metadata$barcode)
MICE$UMAP_1  <- Metadata$UMAP_1
MICE$UMAP_2  <- Metadata$UMAP_2
MICE$cluster <- Metadata$cluster
MICE$sample  <- Metadata$sample
MICE <- join_layers_if_needed(MICE)

counts <- MICE@assays$RNA$counts
rownames(counts) <- toupper(rownames(MICE))
MICE <- CreateSeuratObject(counts = counts, meta.data = MICE@meta.data)

umap <- Metadata[, c("UMAP_1", "UMAP_2")]
rownames(umap) <- Metadata$barcode
MICE[["ref.umap"]] <- CreateDimReducObject(embeddings = as.matrix(umap), key = "UMAP_", assay = DefaultAssay(MICE))

ifnb <- MICE
ifnb <- join_layers_if_needed(ifnb)
ifnb <- SCTransform(ifnb)
ifnb <- RunPCA(ifnb)

set.seed(24022012)
ifnb <- RunHarmony(ifnb, group.by.vars = "sample", dims.use = 1:40)
ifnb <- RunUMAP(ifnb, dims = 1:40, reduction = "harmony", min.dist = 0.3, seed.use = 24022012,
                reduction.name = "umap.harmony", reduction.key = "Uh_", return.model = TRUE)

saveRDS(ifnb, file.path(DERIVED_DIR, "Figure2_F_SAMC_reference.rds"))
cat("Saved: Figure2_F_SAMC_reference.rds\n")
