# ===================== Setup and Arguments =========================
# Read command-line arguments
args = commandArgs(trailingOnly=TRUE)  # Capture command-line arguments
print(args)  # Print the arguments for confirmation

# Set maximum global size for future variables (if using parallel processing)
options(future.globals.maxSize = 100000000 * 1024^2)

# ===================== Set Random Seed =============================

# Set a random seed for reproducibility of results
set.seed(123456789)  # Change this seed to control randomization (use different seeds for variability)

# ===================== Load Libraries ==============================

# Load visualization and data processing libraries
library(viridis)  # Color palettes
library(dplyr)  # Data manipulation
library(Seurat)  # Single-cell RNA-seq analysis
library(readr)  # Reading and writing data
library(cowplot)  # Plotting package
library(patchwork)  # Combining multiple ggplot2 plots
library(ggplot2)  # Plotting package

# Co-expression network analysis libraries
library(igraph)  # Network analysis
library(harmony)  # Harmonization of single-cell RNA-seq data
library(dittoSeq)  # Helper functions for RNA-seq data analysis
library(RColorBrewer)  # Color palettes for plots

# Additional packages for gene annotations and clustering
library(org.Hs.eg.db)  # Gene annotation for human genes
require(tidyverse)  # Includes multiple data manipulation and plotting libraries
library(clusterProfiler)  # Gene enrichment analysis
library(WGCNA)  # Weighted Gene Co-expression Network Analysis
library(hdWGCNA)  # High-dimensional WGCNA
library(scCustomize)  # Customization for Seurat objects

# Libraries for single-cell data manipulation and analysis
library(SCopeLoomR)  # Scalable single-cell data handling
library(dreamlet)  # Differential network analysis for single-cell data
library(SingleCellExperiment)  # Single-cell experiment object
library(scater)  # Single-cell RNA-seq analysis
# library(UCell)  # (commented out) - for single-cell scoring
# library(Nebulosa)  # (commented out) - for density estimation

# Additional packages for plotting and visualization
library(ggplot2)  # Data visualization
library(dplyr)  # Data manipulation
library(RColorBrewer)  # Color palettes for plots
library(cowplot)  # Plotting package
library(ggtree)  # Phylogenetic tree plotting
library(aplot)  # Plotting
library(circlize)  # Circular visualizations
library(ComplexHeatmap)  # Complex heatmaps

# Meta-analysis and statistical analysis packages
library(tidyr)  # Data tidying
library(muscat)  # For single-cell RNA-seq analysis with multi-condition datasets
library(broom)  # Convert statistical results into tidy format
library(tidyverse)  # Data manipulation and visualization
library(metafor)  # Meta-analysis of statistical data

# ===================== Input and Output Setup =====================

# Parse command-line arguments into variables
input_rds <- args[1]  # Path to input RDS file
group_by <- args[2]  # Variable for grouping (e.g., cell types)
only_microglia <- as.logical(args[3])  # Filter for only microglial cells (TRUE/FALSE)
n_of_aggregated <- as.numeric(args[4])  # Number of aggregated samples (e.g., 25)
n_of_aggregated_text <- args[4]  # String version of the number of aggregated samples
threads <- as.numeric(args[5])  # Number of threads for parallel processing
output_pdf <- args[6]  # Path to output PDF file

# ===================== Working Directory Setup =====================

# Define and create output directory
workdir <- "ICH_Stroke/result/7.WGCNA/"  # Directory for WGCNA results
source('utils.R')  # Source utility functions (assumed to be in 'utils.R')
dir.create(workdir)  # Create the working directory if it doesn't exist


#======================= Create Metacells =========================
#======================= Create Metacells =========================
# Read the compressed RDS file containing immune cell data
Immune <- readCRDS('4.Transfer_Labels/Immune_Cells.rds.zstd')
# Subset the 'Immune' Seurat object to include only cells labeled as 'Myeloid' in the 'CellPop' column
Microglia_only <- subset(Immune, subset = CellPop %in% c('Myeloid'))
# Read the predictions CSV file
predictions_csv <- readr::read_csv('/sc/arion/projects/CommonMind/kyriad02/ICH_Stroke/result/Cortex/5.FreshMG_Mapping/MG_predictions_subtype_ADAM.csv')
# Rename the first column to "Cell" for consistency
colnames(predictions_csv)[1] <- "Cell"
# Set row names of the dataframe to the values in the "Cell" column
rownames(predictions_csv) <- predictions_csv$Cell
# Reorder the predictions dataframe according to the column names of the 'Immune' dataframe
predictions_ordered <- predictions_csv[colnames(Immune),]
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
# File path for the UMAP coordinates CSV file
umap_csv_input <- '5.FreshMG_Mapping/MG_predictions_UMAP.csv'
# Read the UMAP coordinates from the CSV file and convert to a data frame
umap <- as.data.frame(readr::read_csv(umap_csv_input))
# Set row names of the data frame to the values in the "CellName" column
rownames(umap) <- umap$CellName
# Remove the first column (assumed to be the "CellName" column)
umap[,1] <- NULL
# Remove the "CellName" column from the data frame
umap$CellName <- NULL
# Reorder the rows of 'umap' to match the column names of the 'Microglia_only' Seurat object
umap <- umap[colnames(Microglia_only),]
# Create a new dimension reduction object ('scanvi_umap') in the 'Microglia_only' Seurat object
# using the UMAP coordinates read from the CSV file
Microglia_only[["scanvi_umap"]] <- CreateDimReducObject(embeddings = as.matrix(umap),  # UMAP coordinates
                                                        key = "UMAP_",  # Key for the new dimension reduction
                                                        assay = DefaultAssay(Microglia_only))  # Assay type



group_by <- 'Subclass'
DefaultAssay(Microglia_only) <- "RNA"
Microglia_only@reductions$pca <- Microglia_only@reductions$pca_regressed_harmony

Microglia_only <- SetupForWGCNA(
    Microglia_only,
    gene_select = "fraction", # the gene selection approach
    fraction = 0.05, # fraction of cells that a gene needs to be expressed in order to be included
    wgcna_name = "hdWGCNA_ICH" # the name of the hdWGCNA experiment
)
# construct metacells  in each group
Microglia_only <- MetacellsByGroups(
    seurat_obj = Microglia_only,
    reduction = 'pca_regressed_harmony',
    group.by = c(group_by, "orig.ident"), # specify the columns in seurat_obj@meta.data to group by
    k = 25, # nearest-neighbors parameter
    ident.group = 'Subclass' # set the Idents of the metacell seurat object
)

# Normalize the gene expression data in the 'Microglia_only' Seurat object
Microglia_only <- NormalizeMetacells(Microglia_only)

# Scale the normalized gene expression data in the 'Microglia_only' Seurat object
# using the variable features identified in the previous step
Microglia_only <- ScaleMetacells(Microglia_only, features = VariableFeatures(Microglia_only))

# Identify variable features in the gene expression data of the 'Microglia_only' Seurat object
Microglia_only <- FindVariableFeatures(object = Microglia_only)

# Perform principal component analysis (PCA) on the variable features of the gene expression data
# in the 'Microglia_only' Seurat object
Microglia_only <- RunPCAMetacells(Microglia_only, features = VariableFeatures(Microglia_only))

set.seed(123456789)
# Perform batch correction using Harmony on the 'Microglia_only' Seurat object
# using the 'orig.ident' as grouping variables
Microglia_only <- RunHarmonyMetacells(Microglia_only, group.by.vars = 'orig.ident',n.seed=123456789)#15678234

# Perform uniform manifold approximation and projection (UMAP) on the 'Microglia_only' Seurat object
# using the Harmony-corrected data and the specified dimensions and parameters
Microglia_only <- RunUMAPMetacells(Microglia_only,
                                    reduction = 'harmony',  # Reduction method
                                    dims = 1:10,            # Dimensions to use
                                    min_dist = 0.5)         # Minimum distance parameter
# Retrieve the metacell object from 'Microglia_only' and perform neighborhood and clustering identification
All_metacell <- GetMetacellObject(Microglia_only) %>%
        FindNeighbors(reduction = "harmony") %>%
        FindClusters(resolution = seq(0.1, 1, 0.1))  # Cluster at different resolutions

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

saveCRDS(All_metacell,paste0(workdir,'Myeloid_Metacells_Subclass_ADAM.rds.ztsd'))

