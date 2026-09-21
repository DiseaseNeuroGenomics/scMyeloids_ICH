# Author: Dimitrios Kyriakis
# ===================== Setup and Arguments =========================
args = commandArgs(trailingOnly=TRUE)
print(args)

options(future.globals.maxSize = 100000000 * 1024^2)
set.seed(123456789)

# ===================== Load Libraries ==============================
library(viridis)
library(dplyr)
library(Seurat)
library(readr)
library(cowplot)
library(patchwork)
library(ggplot2)
library(igraph)
library(harmony)
library(dittoSeq)
library(RColorBrewer)
library(org.Hs.eg.db)
require(tidyverse)
library(clusterProfiler)
library(WGCNA)
library(hdWGCNA)
library(scCustomize)
library(SCopeLoomR)
library(dreamlet)
library(SingleCellExperiment)
library(scater)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(cowplot)
library(ggtree)
library(aplot)
library(circlize)
library(ComplexHeatmap)
library(tidyr)
library(muscat)
library(broom)
library(tidyverse)
library(metafor)

# ===================== Input and Output Setup =====================
input_rds <- args[1]          # Path to input RDS file
group_by <- args[2]           # Variable for grouping (e.g., cell types)
only_microglia <- as.logical(args[3])  # Filter for only microglial cells (TRUE/FALSE)
n_of_aggregated <- as.numeric(args[4]) # Number of aggregated samples (e.g., 25)
n_of_aggregated_text <- args[4]        # String version of the number of aggregated samples
threads <- as.numeric(args[5]) # Number of threads for parallel processing
output_pdf <- args[6]          # Path to output PDF file
predictions_csv <- args[7]     # Path to FreshMG predictions CSV
umap_csv <- args[8]            # Path to UMAP coordinates CSV
output_rds <- args[9]          # Path to output metacell RDS

# ===================== Working Directory Setup =====================
workdir <- dirname(output_rds)
source('utils.R')   # lab helper, not redistributed; see KNOWN_ISSUES.md
dir.create(workdir, recursive = TRUE, showWarnings = FALSE)

#======================= Create Metacells =========================
Immune <- readCRDS(input_rds)
Microglia_only <- subset(Immune, subset = CellPop %in% c('Myeloid'))

# Read the predictions CSV file
pred_df <- readr::read_csv(predictions_csv)
colnames(pred_df)[1] <- "Cell"
rownames(pred_df) <- pred_df$Cell

# Reorder the predictions dataframe according to the column names of the 'Immune' dataframe
predictions_ordered <- pred_df[colnames(Immune),]

# Assign predicted subclass and subtype values to 'Immune' dataframe
Immune$Subclass <- predictions_ordered$predictions_subclass
Immune$Subtype <- predictions_ordered$predictions_subtype

# Replace NA values in 'Subclass' and 'Subtype' with 'CellPop' values
Immune$Subclass[is.na(Immune$Subtype)] <- Immune$CellPop[is.na(Immune$Subtype)]
Immune$Subtype[is.na(Immune$Subtype)] <- Immune$CellPop[is.na(Immune$Subtype)]

# Group similar subtypes under broader subclasses
Immune$Subclass[grep('Adapt',Immune$Subtype)] <- 'Adapt'
Immune$Subclass[grep('ADAM',Immune$Subtype)] <- 'ADAM'
Immune$Subclass[grep('Homeo',Immune$Subtype)] <- 'Homeo'
Immune$Subclass[grep('PVM',Immune$Subtype)] <- 'PVM'
Immune$Subclass[grep('Prolif',Immune$Subtype)] <- 'Prolif'
Immune$Subclass[grep('exAM',Immune$Subtype)] <- 'exAM'

# Read the UMAP coordinates from the CSV file
umap <- as.data.frame(readr::read_csv(umap_csv))
rownames(umap) <- umap$CellName
umap[,1] <- NULL
umap$CellName <- NULL

# Reorder the rows of 'umap' to match the column names of the 'Microglia_only' Seurat object
umap <- umap[colnames(Microglia_only),]

# Create a new dimension reduction object ('scanvi_umap') in the 'Microglia_only' Seurat object
Microglia_only[["scanvi_umap"]] <- CreateDimReducObject(
    embeddings = as.matrix(umap),
    key = "UMAP_",
    assay = DefaultAssay(Microglia_only)
)

group_by <- 'Subclass'
DefaultAssay(Microglia_only) <- "RNA"
Microglia_only@reductions$pca <- Microglia_only@reductions$pca_regressed_harmony

Microglia_only <- SetupForWGCNA(
    Microglia_only,
    gene_select = "fraction",
    fraction = 0.05,
    wgcna_name = "hdWGCNA_ICH"
)

# construct metacells in each group
Microglia_only <- MetacellsByGroups(
    seurat_obj = Microglia_only,
    reduction = 'pca_regressed_harmony',
    group.by = c(group_by, "orig.ident"),
    k = 25,
    ident.group = 'Subclass'
)

Microglia_only <- NormalizeMetacells(Microglia_only)
Microglia_only <- ScaleMetacells(Microglia_only, features = VariableFeatures(Microglia_only))
Microglia_only <- FindVariableFeatures(object = Microglia_only)
Microglia_only <- RunPCAMetacells(Microglia_only, features = VariableFeatures(Microglia_only))

set.seed(123456789)
Microglia_only <- RunHarmonyMetacells(Microglia_only, group.by.vars = 'orig.ident', n.seed=123456789)

Microglia_only <- RunUMAPMetacells(Microglia_only,
                                    reduction = 'harmony',
                                    dims = 1:10,
                                    min_dist = 0.5)

All_metacell <- GetMetacellObject(Microglia_only) %>%
        FindNeighbors(reduction = "harmony") %>%
        FindClusters(resolution = seq(0.1, 1, 0.1))

# Assign metadata from 'Microglia_only' to the 'All_metacell' object based on matching 'orig.ident'
All_metacell$Donor <- All_metacell$orig.ident
All_metacell$mRS <- as.vector(Microglia_only$mRS[match(All_metacell$orig.ident, Microglia_only$orig.ident)])
All_metacell$TSH <- as.vector(Microglia_only$TSH[match(All_metacell$orig.ident, Microglia_only$orig.ident)])
All_metacell$race <- as.vector(Microglia_only$race[match(All_metacell$orig.ident, Microglia_only$orig.ident)])
All_metacell$age <- as.vector(Microglia_only$age[match(All_metacell$orig.ident, Microglia_only$orig.ident)])
All_metacell$sex <- as.vector(Microglia_only$sex[match(All_metacell$orig.ident, Microglia_only$orig.ident)])
All_metacell$mRS_binned <- as.vector(Microglia_only$mRS_binned[match(All_metacell$orig.ident, Microglia_only$orig.ident)])
All_metacell$mRS_3class <- as.vector(Microglia_only$mRS_3class[match(All_metacell$orig.ident, Microglia_only$orig.ident)])
All_metacell$BL <- as.vector(Microglia_only$BL[match(All_metacell$orig.ident, Microglia_only$orig.ident)])
All_metacell$CL <- as.vector(Microglia_only$CL[match(All_metacell$orig.ident, Microglia_only$orig.ident)])
All_metacell$dx <- as.vector(Microglia_only$dx[match(All_metacell$orig.ident, Microglia_only$orig.ident)])
All_metacell$PHV <- as.vector(Microglia_only$PHV[match(All_metacell$orig.ident, Microglia_only$orig.ident)])
All_metacell$PEV <- as.vector(Microglia_only$PEV[match(All_metacell$orig.ident, Microglia_only$orig.ident)])

All_metacell$Clusters <- as.vector(All_metacell$RNA_snn_res.0.4)
All_metacell$Clusters[All_metacell$RNA_snn_res.0.4 %in% c(0,8)] <- 'MTC_1'
All_metacell$Clusters[All_metacell$RNA_snn_res.0.4 %in% c(2)] <- 'MTC_2'
All_metacell$Clusters[All_metacell$RNA_snn_res.0.4 %in% c(5,6)] <- 'MTC_3'
All_metacell$Clusters[All_metacell$RNA_snn_res.0.4 %in% c(3)] <- 'MTC_4'
All_metacell$Clusters[All_metacell$RNA_snn_res.0.4 %in% c(4)] <- 'MTC_5'
All_metacell$Clusters[All_metacell$RNA_snn_res.0.4 %in% c(1)] <- 'MTC_6'
All_metacell$Clusters[All_metacell$RNA_snn_res.0.4==7] <- 'MTC_Prolif'

saveCRDS(All_metacell, output_rds)
message("Saved metacell object to: ", output_rds)
