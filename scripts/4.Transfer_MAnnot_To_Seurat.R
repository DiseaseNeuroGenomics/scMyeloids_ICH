# Set a random seed for reproducibility
set.seed(12345)

# Increase the maximum size of global variables in future package
options(future.globals.maxSize = 100000000 * 1024^2)

# Read in command-line arguments
args = commandArgs(trailingOnly=TRUE)
print(args)

# Load required libraries
library(Seurat)
library(ggplot2)
library(tidyverse)
library(cowplot)
library(gtools)
library(parallel)
library(SeuratDisk)
library(patchwork)

# Define a negation operator for use in filtering
'%notin%' <- Negate('%in%')

# Set a consistent theme for plots
theme_to_add <- theme_bw() + theme(text = element_text(size = 20))

# Source utility functions
source('utils.R')

# ===================== READ ARGUMENTS =========================
# Assign command-line arguments to specific variables
Cortex_data <- args[1]
pca_regressed_harmony <- args[2]
umap_coords <- args[3]
pg_metadata <- args[4]
activation_markers <- args[5]
rds_output <- args[6]
sce_rds_output <- args[7]
output_metadata <- args[8]
# --------------------------------------------------------------

# ================== LOAD Harmonized Projections ================
# Read PCA and UMAP projections from files
PC_sp <- as.data.frame(readr::read_tsv(pca_regressed_harmony))
rownames(PC_sp) <- PC_sp[,1]
PC_sp[,1] <- NULL
PC_sp <- as.matrix(PC_sp)

UMAP_sp <- as.data.frame(readr::read_tsv(umap_coords))
rownames(UMAP_sp) <- UMAP_sp[,1]
UMAP_sp[,1] <- NULL
UMAP_sp <- as.matrix(UMAP_sp)
# --------------------------------------------------------------

# ==================== PREPROCESS DATA =========================
# Load Cortex data and add barcode keys
Cortex_before_Extra_filter <- readRDS(Cortex_data)
Cortex_before_Extra_filter$barcodeKey <- colnames(Cortex_before_Extra_filter)
DefaultAssay(Cortex_before_Extra_filter) <- "RNA"

# Subset Cortex data to match UMAP projections
Cortex <- subset(Cortex_before_Extra_filter, barcodeKey %in% rownames(UMAP_sp))
rm(Cortex_before_Extra_filter)

# Add dimensionality reduction data to Cortex object
Cortex[["pca_regressed_harmony"]] <- CreateDimReducObject(embeddings = PC_sp, key = "pca_regressed_harmony_", assay = DefaultAssay(Cortex))
Cortex[["harmony_umap"]] <- CreateDimReducObject(embeddings = UMAP_sp, key = "harmony_umap_", assay = DefaultAssay(Cortex))

# Run UMAP and clustering analysis
Cortex <- RunUMAP(Cortex, reduction = "pca_regressed_harmony", dims = 1:15, verbose = FALSE, return.model = TRUE, min.dist = 0.5)
Cortex <- FindNeighbors(Cortex, reduction = "pca_regressed_harmony", dims = 1:15, verbose = FALSE)
Cortex <- FindClusters(Cortex, resolution = 0.4)

# Add metadata to Cortex object
pg_meta <- as.data.frame(readr::read_tsv(pg_metadata))
Cortex$predicted_phase <- pg_meta$predicted_phase
Cortex$CC_diff <- pg_meta$CC_diff
Cortex$leiden_labels <- pg_meta$leiden_labels
Cortex$doublet_score <- pg_meta$scrublet_score
Cortex$pred_dbl <- pg_meta$pred_dbl
Cortex$anno <- pg_meta$anno
Cortex$anno_dh <- pg_meta$anno_dh
Cortex$class <- pg_meta$class
Cortex$All.Act.Score <- as.vector(pg_meta["All CNS Cells Act. Score"][,1])
Cortex$Act.Score <- as.vector(pg_meta["Micro/Myeloid Shared Act. Score"][,1])
Cortex$Mic.Ident <- as.vector(pg_meta["Microglial Identity Score"][,1])

# Save UMAP plots to PDF
pdf("1.Cortex_UMAP.pdf")
DimPlot(Cortex, reduction = "umap", group.by = c("race")) 
DimPlot(Cortex, reduction = "umap", group.by = c("sex"))
DimPlot(Cortex, reduction = "umap", group.by = c("donor"))
DimPlot(Cortex, reduction = "umap", group.by = c("class"), label = TRUE)
DimPlot(Cortex, reduction = "umap", group.by = c("anno"), label = TRUE)
DimPlot(Cortex, reduction = "umap", group.by = c("seurat_clusters"), label = TRUE)
FeaturePlot(Cortex, reduction = "umap", "CC_diff")
dev.off()

# Calculate average expression for donors
DefaultAssay(Cortex) <- "RNA"
Cortex_cl <- AverageExpression(Cortex, group.by = "donor", return.seurat = TRUE)
dem_donor <- Cortex_cl@assays$RNA@data
dim(dem_donor)

# ================ Create Activation Score ===================
# Read activation marker genes and process
df <- readr::read_csv(activation_markers)
df$GeneSynbol <- toupper(df$GeneSynbol)
activation_list <- split(df$GeneSynbol, df$Module)

# Add module scores to Cortex object
Cortex <- AddModuleScore(Cortex, activation_list)
Cortex$Act.Score_Allcells <- Cortex$Cluster1
Cortex$Act.Score <- Cortex$Cluster2
Cortex$identity_score <- Cortex$Cluster3
Cortex$Microglial_identity_score <- Cortex$Cluster3

# ================ CELL TYPE ANNOTATION =====================
# Annotate cell types based on clusters
cluster_vec <- as.numeric(Cortex$seurat_clusters) - 1
Cortex$CellType <- cluster_vec
Cortex$CellType[cluster_vec %in% c(2)] <- "MG-Homeo."
Cortex$CellType[cluster_vec %in% c(0)] <- "MG-Active"
Cortex$CellType[cluster_vec %in% c(1)] <- "MG-Inter."
Cortex$CellType[cluster_vec %in% c(3)] <- "Monocytes"
Cortex$CellType[cluster_vec %in% c(4)] <- "Proliferation"
Cortex$CellType[cluster_vec %in% c(5)] <- "T Cells"
Cortex$CellType[cluster_vec %in% c(6)] <- "Astroccytes"
Cortex$CellType[cluster_vec %in% c(7)] <- "Oligodendrocytes"
Cortex$CellType[cluster_vec %in% c(8)] <- "Murel"
Cortex$CellType[cluster_vec %in% c(9)] <- "B Cells"
Cortex$CellType[cluster_vec %in% c(10)] <- "exclude"

# Simplify cell populations
Cortex$CellPop <- Cortex$CellType
Cortex$CellPop[Cortex$CellType %in% c("MG-Homeo.", "MG-Inter.", "MG-Active")] <- "Microglia"
Cortex$CellPop[Cortex$class == "Monocytes"] <- "Monocytes"
Cortex <- subset(Cortex, subset = CellType != 'exclude')

# Save cell cycle gene plots
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
pdf("Cell_Cycle_Genes.pdf", width = 12)
FeaturePlot(Cortex, features = c(s.genes, g2m.genes), order = TRUE, raster = TRUE)
VlnPlot(Cortex, features = c(s.genes, g2m.genes), sort = TRUE, flip = TRUE, stack = TRUE)
dev.off()

# Save activation score scatter plots
pdf("1.Activation_Score_IN_cortex.pdf", width = 12)
FeatureScatter(Cortex, "Act.Score_Allcells", "identity_score", shape.by = "class", group.by = "class", raster = TRUE, pt.size = 0.3) +
  xlab("Activation score") + ylab("Microglial identity score") +
  geom_hline(yintercept = -0.15, linetype = 'dashed', col = 'black', size = 0.7) +
  geom_vline(xintercept = 0.25, linetype = 'dashed', col = 'black', size = 0.7) +
  theme_bw() + theme(text = element_text(size = 20))
FeatureScatter(Cortex, "Act.Score_Allcells", "identity_score", shape.by = "CellType", group.by = "CellType", raster = TRUE, pt.size = 0.3) +
  xlab("Activation score") + ylab("Microglial identity score") +
  geom_hline(yintercept = -0.15, linetype = 'dashed', col = 'black', size = 0.7) +
  geom_vline(xintercept = 0.25, linetype = 'dashed', col = 'black', size = 0.7) +
  theme_bw() + theme(text = element_text(size = 20))
dev.off()

# Save combined scatter and feature plots
pdf("2.Scatter_Activation_Score_IN_cortex.pdf", width = 13, height = 10)
p1 <- DimPlot(Cortex, group.by = "CellType", reduction = "umap", label = TRUE) + theme_to_add
p3 <- FeaturePlot(Cortex, reduction = "umap", features = c("Act.Score_Allcells"), order = TRUE) + theme_to_add + ggtitle("Activation score")
p4 <- FeaturePlot(Cortex, reduction = "umap", features = c("identity_score"), order = TRUE) + theme_to_add + ggtitle("Microglial identity score")
p2 <- FeatureScatter(Cortex, "Act.Score_Allcells", "identity_score", group.by = "CellType", raster = TRUE) +
  xlab("Activation score") + ylab("Microglial identity score") +
  geom_hline(yintercept = -0.15, linetype = 'dashed', col = 'black', size = 0.7) +
  geom_vline(xintercept = 0.25, linetype = 'dashed', col = 'black', size = 0.7) +
  theme_bw() + theme(text = element_text(size = 20))
p1 + p2 + p3 + p4 + plot_layout(ncol = 2, guides = "collect")
dev.off()

# Save Cortex object and metadata
saveCRDS(Cortex, rds_output)
write_tsv(as.data.frame(Cortex@meta.data), output_metadata)

# Convert Cortex object to SingleCellExperiment and save
sce_r <- as.SingleCellExperiment(Cortex)
saveCRDS(sce_r, sce_rds_output)
# ----------------------------------------------------------------------
