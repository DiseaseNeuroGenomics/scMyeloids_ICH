# Author: Dimitrios Kyriakis
args = commandArgs(trailingOnly=TRUE)
print(args)
options(future.globals.maxSize = 100000000 * 1024^2)
set.seed(123456789)
print(args)

# ========================= Libraries ===================================
library(dplyr)
library(tidyr)
library(tidyverse)
library(readr)
library(Seurat)
library(SingleCellExperiment)
library(scater)
library(scuttle)
library(Seurat)
library(igraph)
library(harmony)
library(WGCNA)
library(hdWGCNA)
library(clusterProfiler)
library(SCopeLoomR)
library(dreamlet)
library(ggplot2)
library(viridis)
library(cowplot)
library(patchwork)
library(RColorBrewer)
library(ComplexHeatmap)
library(ggtree)
library(aplot)
library(circlize)
library(org.Hs.eg.db)
library(muscat)
library(broom)
library(metafor)

'%notin%' <- Negate('%in%')

# ========================= Input arguments ===================================
input_rds <- args[1]                # Input RDS file
group_by <- args[2]                 # Grouping variable
only_microglia <- as.logical(args[3])  # Flag for selecting microglia data
n_of_aggregated <- as.numeric(args[4]) # Number of aggregated data points
n_of_aggregated_text <- args[4]         # Text version for output
threads <- as.numeric(args[5])      # Number of threads for parallel processing
output_pdf <- args[6]               # Output PDF filename
output_tsv <- args[7]               # Output TSV for DEG table

# ========================= Working Directory ===================================
workdir <- dirname(input_rds)
source('scripts/utils.R')
dir.create(workdir, recursive = TRUE, showWarnings = FALSE)

# Load the metacell data from a compressed RDS file
Metacells <- readCRDS(input_rds)

# Create a new cluster grouping with three major categories based on original clusters
Metacells$ThreeClusters <- as.vector(Metacells$Clusters)
Metacells$ThreeClusters[Metacells$Clusters %in% c("MTC_3")] <- "Homeostatic_MTC"
Metacells$ThreeClusters[Metacells$Clusters %in% c("MTC_1", "MTC_2")] <- "Lineage2"
Metacells$ThreeClusters[Metacells$Clusters %in% c("MTC_4", "MTC_5", "MTC_6", "MTC_7", "MTC_Prolif")] <- "Lineage1"

# Create a new cluster grouping with two major categories based on original clusters
Metacells$TwoClusters <- as.vector(Metacells$Clusters)
Metacells$TwoClusters[Metacells$Clusters %in% c("MTC_1", "MTC_2", "MTC_3")] <- "MTC_123"
Metacells$TwoClusters[Metacells$Clusters %in% c("MTC_4", "MTC_5", "MTC_6")] <- "MTC_456"

# Subset data to include only cells from the two main cluster groups ("MTC_123" and "MTC_456")
subset_two_cl <- subset(Metacells, subset = TwoClusters %in% c("MTC_123", "MTC_456"))

# Set cluster identities to "TwoClusters" for further analysis
Idents(subset_two_cl) <- "TwoClusters"

# Identify DEGs between the two clusters with no log fold-change threshold
DEGs1 <- FindAllMarkers(subset_two_cl, only.pos = T, logfc.threshold = 0, return.thresh = 1)

# Calculate fold-change between "MTC_456" and "MTC_123" groups with minimal filtering
fold_ch_res <- FoldChange(
    subset_two_cl,
    ident.1 = 'MTC_456',
    ident.2 = 'MTC_123',
    logfc.threshold = 0,
    slot = "data",
    min.pct = 0.1,
    min.cells.group = 3,
    pseudocount.use = 1,
    mean.fxn = NULL,
    fc.name = NULL,
    base = 2,
    return.thresh = 1
)

# Calculate standard deviation of expression values for cells in "MTC_456" cluster
obj <- subset_two_cl@assays$RNA@data[, subset_two_cl$TwoClusters == 'MTC_456']
row_sd <- apply(obj, 1, sd)
fold_ch_res$sd.1 <- row_sd

# Calculate standard deviation of expression values for cells in "MTC_123" cluster
obj <- subset_two_cl@assays$RNA@data[, subset_two_cl$TwoClusters == 'MTC_123']
row_sd <- apply(obj, 1, sd)
fold_ch_res$sd.2 <- row_sd

# Calculate average expression levels for each gene across the "TwoClusters" groups
averageMT <- as.data.frame(AverageExpression(subset_two_cl, group.by = 'TwoClusters'))

# Merge DEG data with fold-change results by row names
merged_df <- merge(DEGs1, fold_ch_res, by = "row.names", all = TRUE)
head(merged_df)
rownames(merged_df) <- merged_df$Row.names
merged_df[, 1] <- NULL

# Merge the merged DEG and fold-change data with average expression data
merged_df2 <- merge(merged_df, averageMT, by = "row.names", all = TRUE)
rownames(merged_df2) <- merged_df2$Row.names
merged_df2[, 1] <- NULL

# Calculate standard error for fold-change values based on sample sizes
se <- sqrt((merged_df2$sd.1^2 / 13313) + (merged_df2$sd.2^2 / 15274))
merged_df2$se <- se

# Compute t-statistics based on log2 fold-change values and standard errors
merged_df2$log2FC <- ifelse(merged_df2$cluster == 'MTC_456', merged_df2$avg_log2FC.x, -1 * merged_df2$avg_log2FC.x)
merged_df2$tstat <- merged_df2$log2FC / merged_df2$se

# Remove rows with missing values and arrange genes by descending log2 fold-change
clean_DEGs <- merged_df2 %>% na.omit()
clean_DEGs$gene <- rownames(clean_DEGs)
clean_DEGs <- clean_DEGs %>% arrange(desc(clean_DEGs$log2FC))

# Save the final table of DEGs to a TSV file
readr::write_tsv(clean_DEGs, output_tsv)
message("Saved DEG table to: ", output_tsv)
