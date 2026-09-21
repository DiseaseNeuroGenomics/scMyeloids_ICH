# ============================================================================
# Figure2_G_00_build_ExtendedGenes_loom.R
# Project: Decoding Immune Landscapes and Therapeutic Targets in Intracerebral Hemorrhage
# Author:  Dimitrios Kyriakis
# ============================================================================

# ==============================================================================
# Figure 2, Panel G (input build) — rebuild the "ExtendedGenes" input
# expression loom that feeds the SCENIC GRN pipeline, reproducing
# 4.Figure2_G_Prepare_For_SCENIC.ipynb's gene-curation logic against our own
# cached Metacells object.
#
# Distinct output name (2026_08_ExtGenes_input.loom) from the original
# 2024_04_01_Myeloid_Metacells_ExtendedGenes.loom on purpose — that file is
# being kept as a comparison baseline, not overwritten.
#
# Gene set = union of:
#   - genes with a resolvable ENTREZID/ENSEMBL/GENENAME (bitr on the PCA
#     feature-loading genes)
#   - genes present in the cisTarget v10 motif-ranking database
#   - genes present in the motif annotation table
#   - all known human TFs (allTFs_hg38.txt) that are present in the dataset
#   - top 5000 variable features (FindVariableFeatures)
# matching the original notebook's cells 7-19 exactly; the notebook's own
# `data <- read_feather(feather_file)` in its cell 10 is dead code (read but
# never used for anything downstream there either) and is skipped here.
#
# Input:  DERIVED_DIR/Figure2_Metacells.rds       (Figure2_A_umap_metacells.R)
#         database/allTFs_hg38.txt
#         database/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl
# Output: DERIVED_DIR/2026_08_ExtGenes_input.loom
# ==============================================================================

source("../../00_setup.R")
library(SCopeLoomR)
library(org.Hs.eg.db)
library(clusterProfiler)

# add_cell_annotation(): inlined from aertslab/SCENIC (not SCopeLoomR, which is
# why the currently-installed SCopeLoomR-only env doesn't have it) — a thin
# wrapper around SCopeLoomR's own add_col_attr()/get_cell_ids(), both already
# available, so not worth installing the full SCENIC package for. Verbatim
# from SCENIC's source, confirmed by the user directly.
add_cell_annotation <- function(loom, cellAnnotation) {
  cellAnnotation <- data.frame(cellAnnotation)
  if (any(c("nGene", "nUMI") %in% colnames(cellAnnotation))) {
    warning("Columns 'nGene' and 'nUMI' will not be added as annotations to the loom file.")
    cellAnnotation <- cellAnnotation[, colnames(cellAnnotation) != "nGene", drop = FALSE]
    cellAnnotation <- cellAnnotation[, colnames(cellAnnotation) != "nUMI", drop = FALSE]
  }

  if (ncol(cellAnnotation) <= 0) stop("The cell annotation contains no columns")
  if (!all(get_cell_ids(loom) %in% rownames(cellAnnotation))) stop("Cell IDs are missing in the annotation")

  cellAnnotation <- cellAnnotation[get_cell_ids(loom), , drop = FALSE]
  for (cn in colnames(cellAnnotation)) {
    add_col_attr(loom = loom, key = cn, value = cellAnnotation[, cn])
  }

  invisible(loom)
}

DATABASE_DIR <- file.path(BASE_DIR, "database")
output_loom  <- file.path(DERIVED_DIR, "2026_08_ExtGenes_input.loom")

Metacells <- readRDS(file.path(DERIVED_DIR, "Figure2_Metacells.rds"))

# ---- Genes with a resolvable Entrez/Ensembl ID, from the PCA feature loadings ----
var_genes <- rownames(Metacells@reductions$pca@feature.loadings)
idsFromSYMBOL <- bitr(var_genes,
                       fromType = "SYMBOL", toType = c("ENTREZID", "ENSEMBL", "GENENAME"),
                       OrgDb = org.Hs.eg.db, drop = FALSE)
idsFromSYMBOL_full <- na.omit(idsFromSYMBOL)
idsFromSYMBOL_full_unique <- idsFromSYMBOL_full[!duplicated(idsFromSYMBOL_full$SYMBOL), ]
rownames(idsFromSYMBOL_full_unique) <- idsFromSYMBOL_full_unique$SYMBOL

# ---- TFs and motif-database genes present in this dataset ----
allTFs_path <- file.path(DATABASE_DIR, "allTFs_hg38.txt")
motifs_path <- file.path(DATABASE_DIR, "motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl")

AllTFs <- readr::read_tsv(allTFs_path, col_names = FALSE)
TF_vector <- unique(AllTFs$X1)
Tfs_in_dataset <- TF_vector[TF_vector %in% rownames(Metacells)]

motifs <- readr::read_tsv(motifs_path)
motifs_genes <- unique(motifs$gene_name)
motifs_in_dataset <- motifs_genes[motifs_genes %in% rownames(Metacells)]

cat("TFs in dataset:", length(Tfs_in_dataset), "\n")
cat("Motif-database genes in dataset:", length(motifs_in_dataset), "\n")

# ---- Top 5000 variable features ----
Metacells <- FindVariableFeatures(Metacells, nfeatures = 5000)
features <- VariableFeatures(Metacells)

# ---- Union of all gene sources ----
features_ext <- unique(c(idsFromSYMBOL_full_unique$SYMBOL, features))
genes_to_add <- unique(c(motifs_in_dataset, Tfs_in_dataset))
genes_to_loom <- unique(c(genes_to_add, features_ext))
cat("Total genes in ExtendedGenes loom:", length(genes_to_loom), "\n")

# ---- Build and write the loom ----
exprMat_filtered <- Metacells@assays$RNA@data[genes_to_loom, ]
cellInfo <- Metacells@meta.data

loom <- build_loom(output_loom, dgem = exprMat_filtered)
loom <- add_cell_annotation(loom, cellInfo)
close_loom(loom)

cat("Saved:", output_loom, "\n")
