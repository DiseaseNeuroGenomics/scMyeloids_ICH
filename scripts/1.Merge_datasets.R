# Author: Dimitrios Kyriakis
# Set options and seed for reproducibility
options(future.globals.maxSize = 1000000 * 1024^2)  # Increase memory limit for large data processing
set.seed(12345)  # Set random seed for reproducibility

# ============== Read Arguments ===================
args <- commandArgs(trailingOnly = TRUE)  # Read command-line arguments
print(args)

# Check if at least one argument is provided, else stop with an error message
if (length(args) == 0) {
    stop("At least one argument must be supplied (input file).", call. = FALSE)
} else {
    print("Arguments Passed")  # Indicate that arguments are passed successfully
}

# Parse input arguments
n_args <- length(args)
print(n_args)
index_rds <- n_args - 5  # Determine the index for separating RDS files from metadata and other inputs

# Extract input paths and parameters
rds_files_path <- args[1:index_rds]  # Paths to RDS files
print(rds_files_path)
metada_file <- args[index_rds + 1]  # Metadata file path
mc_cores <- as.integer(args[index_rds + 2])  # Number of cores for parallel processing
hdad_output <- args[index_rds + 3]  # Output directory for HDAD results
rds_output <- args[index_rds + 4]  # Output file for merged RDS object
merged_qc <- args[n_args]  # QC plot output path

# ================ Load Libraries ==============
library(Seurat)
library(SeuratDisk)
library(harmony)
library(parallel)
library(tidyverse)
library(ggplot2)
library(cowplot)
library(gtools)
library(patchwork)
# ----------------------------------------------

# List of Seurat objects to merge
seurat_list <- rds_files_path
print(seurat_list)

# Load metadata and ensure proper formatting
metadata <- readr::read_tsv(metada_file)  # Read metadata from the provided file
metadata <- as.data.frame(metadata %>% dplyr::rename_all(make.names))  # Standardize column names
print(colnames(metadata))  # Print column names to verify proper loading

# Load RDS files into a list of Seurat objects
seurat_list_obj <- mclapply(
    X = seurat_list,
    FUN = function(seur) { readRDS(seur) },
    mc.cores = mc_cores
)

# Merge Seurat objects into a single dataset
Cortex <- merge(
    seurat_list_obj[[1]],
    y = seurat_list_obj[2:length(seurat_list_obj)],
    project = "ICH_Cortex"
)

# Save the merged Seurat object
saveRDS(Cortex, rds_output)

# Save the merged object in H5Seurat format and convert to AnnData (h5ad) format
SaveH5Seurat(Cortex, filename = file.path(dirname(rds_output), "Merged_Cortex.h5Seurat"), overwrite = TRUE)
Convert(file.path(dirname(rds_output), "Merged_Cortex.h5Seurat"), dest = "h5ad")
