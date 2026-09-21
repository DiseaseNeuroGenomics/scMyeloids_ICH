# Author: Dimitrios Kyriakis
# Set a seed for reproducibility of results
set.seed(1234)

# Increase the maximum size for global variables to handle large datasets
options(future.globals.maxSize = 1000000 * 1024^2)

# Retrieve command-line arguments passed to the script
args = commandArgs(trailingOnly = TRUE)

# Load required libraries for data manipulation, visualization, and analysis
library(Seurat)         # For single-cell RNA sequencing analysis
library(ggplot2)        # For data visualization
library(sctransform)    # For single-cell normalization and variance stabilization
library(tidyverse)      # For data manipulation and visualization
library(cowplot)        # For creating compound plots
library(dropestr)       # For processing drop-seq data
library(dittoSeq)       # For plotting and visualization of single-cell data
library(pbapply)        # For progress-bar-enabled apply functions
library(parallel)       # For parallel computation

# Define a custom negation operator for `%in%`
`%notin%` <- Negate(`%in%`)

# Source external functions for preprocessing
source("0.Functions.R")

# Assign command-line arguments to variables for further use
h5_file <- args[1]                 # Path to the HDF5 file containing gene expression data
metada_file <- args[2]             # Path to the metadata file
origin_workdir <- args[3]          # Path to the working directory
sample <- args[4]                  # Sample name
min.molecules.per.gene <- as.integer(args[5])  # Minimum number of molecules required per gene
min.cells <- 10                    # Minimum number of cells required for filtering
remove_mitochondrial <- args[6]    # Whether to remove mitochondrial genes
remove_ribosomal <- args[7]        # Whether to remove ribosomal genes
remove_mitocard <- args[8]         # Whether to remove mitochondrial CARD genes
remove_noncoding <- args[9]        # Whether to remove noncoding genes
rdsoutout <- args[10]              # Output path for the RDS file
reported_MAD_fig <- args[11]       # Path for the output MAD figure
reported_scrublet_fig <- args[12]  # Path for the output Scrublet figure
reported_After_MAD_fig <- args[13] # Path for the post-MAD filtering figure

# Load the metadata file and format it
medical_metadata <- readr::read_tsv(metada_file)  # Read metadata as a tab-separated file
medical_metadata <- as.data.frame(medical_metadata %>% dplyr::rename_all(make.names)) # Ensure column names are syntactically valid in R

# Check if the output RDS file already exists
if (file.exists(rdsoutout)) {
    # If the RDS file exists, print a message and skip processing
    print("Prepared")
} else {
    # Create necessary directories for the analysis
    dir.create('result')                                      # Create a result directory
    dir.create(paste0(origin_workdir))                        # Create the base working directory
    dir.create(paste0(origin_workdir, '0.Precheck/'))         # Create a precheck directory for the sample
    workdir <- paste0(origin_workdir, '0.Precheck/', sample, '/')  # Set up the specific work directory for the sample
    dir.create(workdir)

    # Print the sample name for logging purposes
    print(sample)

    # Read the 10X Genomics data using a custom function
    DEM <- read_10x_fun(h5_file = h5_file)

    # Create a Seurat object for single-cell RNA-seq analysis
    seur <- CreateSeuratObject(
        counts = DEM,
        project = sample,
        min.cells = 10,
        min.features = min.molecules.per.gene
    )

    # Add a new metadata column to the Seurat object for the percentage of mitochondrial genes
    seur[["percent.mt"]] <- PercentageFeatureSet(seur, pattern = "^MT-")

    # Extract raw RNA counts from the Seurat object
    DEM <- seur@assays$RNA@counts

    # =================================== GENES =====================================================
    # Identify mitochondrial and ribosomal genes in the dataset using regex
    mt_pos <- grep("^Mt-|^RPL|^RPS", rownames(DEM), ignore.case = TRUE)

    # Load the MitoCarta gene list, which contains mitochondrial genes
    mitocard_file <- readr::read_tsv(args[14])  # Pass MitoCarta file as argument
    mitocard <- unlist(
        lapply(as.vector(mitocard_file$Genes), function(x) {
            as.list(strsplit(x, ", ")[[1]])
        })
    ) %>% na.omit()  # Flatten and remove NA values

    # Load coding and non-coding gene information
    coding_genes <- readr::read_tsv(args[15], col_names = FALSE)  # Pass gencode file as argument
    non_coding_genes <- coding_genes$X2[coding_genes$X3 != "protein_coding"]  # Extract non-coding genes

    # Calculate the percentage of mitochondrial and ribosomal genes for each cell
    seur[["percent.mt"]] <- PercentageFeatureSet(seur, pattern = "^MT-")  # Mitochondrial genes
    seur[["percent.ribo"]] <- PercentageFeatureSet(seur, pattern = "^RPL|^RPS")  # Ribosomal genes

    # Store the percentage metrics for later use
    percent.mt <- seur$percent.mt
    percent.ribo <- seur$percent.ribo

    # Print the dimensions of the raw data matrix
    print(dim(DEM))

    # Optionally remove mitochondrial genes from the dataset
    if (remove_mitochondrial) {
        print("Remove Mitochondrial Genes")
        mitochondrial_pos <- grep("^Mt-", rownames(DEM), ignore.case = TRUE)  # Identify mitochondrial genes
        DEM <- DEM[-mitochondrial_pos, ]  # Remove mitochondrial genes
        print(paste0("Number of Mitochondrial genes removed: ", length(mitochondrial_pos)))
    }

    # Optionally remove ribosomal genes from the dataset
    if (remove_ribosomal) {
        print("Remove Ribosomal Genes")
        ribosomal_pos <- grep("^RPL|^RPS", rownames(DEM), ignore.case = TRUE)  # Identify ribosomal genes
        DEM <- DEM[-ribosomal_pos, ]  # Remove ribosomal genes
        print(paste0("Number of Ribosomal genes removed: ", length(ribosomal_pos)))
    }

    # Optionally remove non-coding genes from the dataset
    if (remove_noncoding) {
        print("Remove Noncoding Genes")
        noncoding_pos <- match(non_coding_genes, rownames(DEM)) %>% na.omit()  # Match non-coding genes
        DEM <- DEM[-noncoding_pos, ]  # Remove non-coding genes
        print(paste0("Number of Noncoding genes removed: ", length(noncoding_pos)))
    }

    # Optionally remove MitoCarta genes from the dataset
    if (remove_mitocard) {
        print("Remove Mitocard Genes")
        mitocard_pos <- match(mitocard, rownames(DEM)) %>% na.omit()  # Match Mitocard genes
        DEM <- DEM[-mitocard_pos, ]  # Remove Mitocard genes
        print(paste0("Number of Mitocard genes removed: ", length(mitocard_pos)))
    }

    # Combine all removed gene positions into a unique list for reporting
    genes_to_remove <- unique(c(mitochondrial_pos, ribosomal_pos, mitocard_pos, noncoding_pos))
    length(genes_to_remove)  # Report the total number of unique genes removed

    # Store the filtered gene expression matrix after removing specified genes
    DEM_pass_genes <- DEM

    # Print the dimensions of the filtered gene expression matrix
    print(dim(DEM_pass_genes))

    # ============================ CELLS ==============================================
    # Create a new Seurat object using the filtered gene expression matrix
    seur <- CreateSeuratObject(
        counts = DEM_pass_genes,
        project = sample,
        min.cells = 0,
        min.features = 0
    )

    # Add metadata for mitochondrial and ribosomal gene percentages to the Seurat object
    seur$percent.mt <- percent.mt  # Percentage of mitochondrial genes
    seur$percent.ribo <- percent.ribo  # Percentage of ribosomal genes

    # Extract metadata and key metrics from the Seurat object
    metadata <- seur@meta.data  # Metadata of the Seurat object
    n_umi <- seur$nCount_RNA    # Number of unique molecular identifiers (UMIs) per cell
    n_genes <- seur$nFeature_RNA  # Number of expressed genes per cell
    n_pc.mit <- seur$percent.mt  # Percentage of mitochondrial genes per cell

    # Calculate MAD-based cutoff limits for the number of UMIs and genes
    n_umi_lims <- calc_mad_cutoff(num_vec = n_umi)
    n_genes_lims <- calc_mad_cutoff(num_vec = n_genes)

    # Generate a PDF report visualizing the MAD cutoffs
    pdf(reported_MAD_fig)

    # Plot histogram of library sizes (log10-scaled) with MAD cutoff lines
    hist(log10(n_umi), col = "lightblue", xlim = c(2, 5))
    abline(v = log10(n_umi_lims[1]), col = "red", lwd = 3, lty = 2)  # Lower limit
    abline(v = log10(n_umi_lims[2]), col = "red", lwd = 3, lty = 2)  # Upper limit

    # Plot histogram of the number of genes (log10-scaled) with MAD cutoff lines
    hist(log10(n_genes), col = "lightblue", xlim = c(2, 4))
    abline(v = log10(n_genes_lims[1]), col = "red", lwd = 3, lty = 2)  # Lower limit
    abline(v = log10(n_genes_lims[2]), col = "red", lwd = 3, lty = 2)  # Upper limit

    # Scatter plot of library size vs. number of expressed genes with cutoff lines
    ggplot(seur@meta.data, aes(x = nCount_RNA, y = nFeature_RNA)) +
        geom_point() +
        geom_vline(xintercept = n_umi_lims, color = "red") +  # Vertical MAD cutoffs for UMIs
        geom_hline(yintercept = n_genes_lims, color = "red") +  # Horizontal MAD cutoffs for genes
        theme_bw()

    dev.off()  # Close the PDF device

    # Identify cells to remove based on UMI, gene, and mitochondrial percentage thresholds
    remove_umi <- names(n_umi)[n_umi < n_umi_lims[1] | n_umi > n_umi_lims[2]]  # Cells with outlier UMIs
    remove_genes <- names(n_genes)[n_genes < n_genes_lims[1] | n_genes > n_genes_lims[2]]  # Cells with outlier genes
    remove_mit <- colnames(seur)[n_pc.mit > 10]  # Cells with >10% mitochondrial gene content

    # Combine all cells flagged for removal into a unique list
    cells_to_remove <- unique(c(remove_umi, remove_genes, remove_mit))

    # Filter out flagged cells from the gene expression matrix
    DEM_pass_mad <- DEM_pass_genes[, colnames(DEM_pass_genes) %notin% cells_to_remove]

    # Update metadata to exclude removed cells
    metadata <- metadata[colnames(seur) %notin% cells_to_remove, ]

    # Print dimensions of the filtered datasets
    dim(DEM_pass_genes)  # Dimensions of the initial filtered gene matrix
    dim(DEM_pass_mad)    # Dimensions of the gene matrix after MAD-based filtering

    # Print the total number of cells removed
    length(cells_to_remove)

    # =================================== scrublet =====================================================
    scrublet_doublet_info <- GetScrubletScores(mat=as.matrix(DEM_pass_mad),sample_name=sample,min.molecules.per.gene=min.molecules.per.gene , method="scrublet",workdir=workdir)

    # =================================== CreateSeuratObject =====================================================
    # Add scrublet doublet information to the metadata
    metadata$is_doublet <- scrublet_doublet_info$is_doublet  # Flag indicating whether a cell is a doublet
    metadata$scrublet_score <- scrublet_doublet_info$score  # Scrublet-generated score for doublet likelihood

    # Ensure that the rownames of metadata match the column names of the filtered gene expression matrix
    table(rownames(metadata) == colnames(DEM_pass_mad))

    # Print the dimensions of metadata and filtered expression matrix
    dim(metadata)  # Dimensions of metadata
    dim(DEM_pass_mad)  # Dimensions of the filtered expression matrix

    # Create a new Seurat object with the filtered gene expression matrix and updated metadata
    seur <- CreateSeuratObject(
        counts = DEM_pass_mad,
        project = sample,
        min.cells = 0,
        min.features = 0,
        meta.data = metadata
    )

    # Generate violin plots for key metrics grouped by sample identity
    p1 <- VlnPlot(
        seur,
        features = c("nCount_RNA", "nFeature_RNA", "percent.mt", "percent.ribo"),
        group.by = "orig.ident"
    )

    # Normalize the data using log normalization
    seur <- NormalizeData(
        seur,
        normalization.method = "LogNormalize",
        scale.factor = 10000
    )

    # Identify highly variable features for downstream analysis
    seur <- FindVariableFeatures(
        seur,
        selection.method = "vst",
        nfeatures = 2000
    )

    # Retrieve all gene names from the Seurat object
    all.genes <- rownames(seur)

    # Scale the data to ensure uniform variance across features
    seur <- ScaleData(seur, verbose = FALSE)

    # ====== Cell Cycle Scoring ======
    s.genes <- cc.genes$s.genes
    g2m.genes <- cc.genes$g2m.genes

    # Perform PCA (Principal Component Analysis) using variable features
    seur <- RunPCA(seur, features = VariableFeatures(object = seur))

    # Identify nearest neighbors and clusters
    seur <- FindNeighbors(seur, dims = 1:20)  # Find neighbors using the top 20 principal components
    seur <- FindClusters(seur)  # Identify clusters

    # Run UMAP for dimensionality reduction and visualization
    seur <- RunUMAP(seur, dims = 1:20)

    # Generate a PDF for UMAP and other visualizations
    pdf(paste0(workdir, sample, "_UMAP.pdf"))

    # Plot key visualizations
    plot(VlnPlot(seur, features = c("nCount_RNA", "nFeature_RNA", "percent.mt", "percent.ribo"), group.by = "orig.ident"))  # Violin plots
    plot(DimPlot(seur, group.by = "is_doublet"))  # UMAP plot colored by doublet status
    plot(DimPlot(seur, label = TRUE, label.box = TRUE))  # UMAP plot with cluster labels
    if (sample != "CK-049-buffy") {
        plot(DimPlot(seur, group.by = "Phase"))  # UMAP plot colored by cell cycle phase
    }
    plot(DimPlot(seur, reduction = "pca", label = TRUE, label.box = TRUE))  # PCA plot with cluster labels
    plot(FeaturePlot(seur, features = c("nCount_RNA", "nFeature_RNA", "percent.mt", "percent.ribo"), order = TRUE))  # Feature plots for key metrics
    plot(FeaturePlot(seur, features = c("IL10", "CD36", "TGFB1", "PPARG"), order = TRUE))  # Feature plots for specific genes
    plot(FeaturePlot(seur, features = c("IL1B", "IL16", "CD68", "IFNG"), order = TRUE))  # Feature plots for additional genes

    dev.off()  # Close the PDF device

    # Extract UMI counts and gene counts for MAD-based cutoff calculation
    n_umi <- seur$nCount_RNA  # UMI counts per cell
    n_genes <- seur$nFeature_RNA  # Gene counts per cell

    # Calculate MAD-based cutoff limits for UMIs and gene counts
    n_umi_lims <- calc_mad_cutoff(num_vec = n_umi)
    n_genes_lims <- calc_mad_cutoff(num_vec = n_genes)

    # Generate a PDF report for MAD-based filtering results
    pdf(reported_After_MAD_fig)

    # Plot histogram of library sizes (log10-scaled) with MAD cutoff lines
    hist(log10(n_umi), col = "lightblue", xlim = c(2, 5))
    abline(v = log10(n_umi_lims[1]), col = "red", lwd = 3, lty = 2)  # Lower limit
    abline(v = log10(n_umi_lims[2]), col = "red", lwd = 3, lty = 2)  # Upper limit

    # Plot histogram of the number of genes (log10-scaled) with MAD cutoff lines
    hist(log10(n_genes), col = "lightblue", xlim = c(2, 4))
    abline(v = log10(n_genes_lims[1]), col = "red", lwd = 3, lty = 2)  # Lower limit
    abline(v = log10(n_genes_lims[2]), col = "red", lwd = 3, lty = 2)  # Upper limit

    dev.off()  # Close the PDF device


    # ============================== TRANSFER METADATA =================================================
    seur$donor <- seur$orig.ident
    seur$race <- medical_metadata$race[match(seur$donor,medical_metadata$donor)]
    seur$age <- medical_metadata$age[match(seur$donor,medical_metadata$donor)]
    seur$etiology <- medical_metadata$etiology[match(seur$donor,medical_metadata$donor)]
    seur$Viability <- medical_metadata$Viability[match(seur$donor,medical_metadata$donor)]
    seur$sex <- medical_metadata$sex[match(seur$donor,medical_metadata$donor)]
    seur$pathology <- medical_metadata$pathology[match(seur$donor,medical_metadata$donor)]
    seur$dx <- medical_metadata$dx[match(seur$donor,medical_metadata$donor)]
    seur$DateOfProcedure <- medical_metadata$DateOfProcedure[match(seur$donor,medical_metadata$donor)]
    seur$PMI_mins <- medical_metadata$PMI_mins[match(seur$donor,medical_metadata$donor)]
    seur$date_processed <- medical_metadata$date_processed[match(seur$donor,medical_metadata$donor)]
    seur$DfCtB <- medical_metadata$Distance.from.the.Clot.to.the.Biopsy..mm[match(seur$donor,medical_metadata$donor)]
    seur$NIHSS <- medical_metadata$NIHSS.on.admission[match(seur$donor,medical_metadata$donor)]
    seur$PMH <- medical_metadata$Past.Medical.History[match(seur$donor,medical_metadata$donor)]
    seur$PHV <- medical_metadata$Preoperative.Hematoma.Volume..mLs.[match(seur$donor,medical_metadata$donor)]
    seur$PEV <- medical_metadata$Preoperative.Edema.Volume..mLs.[match(seur$donor,medical_metadata$donor)]
    seur$BL <- medical_metadata$Biopsy.Location[match(seur$donor,medical_metadata$donor)]
    seur$CL <- medical_metadata$Clot.Location[match(seur$donor,medical_metadata$donor)]
    seur$TSH <- medical_metadata$Time.Since.Hemorrhage..hrs.[match(seur$donor,medical_metadata$donor)]
    seur$mRS <-medical_metadata$mRS.at.6.Months[match(seur$donor,medical_metadata$donor)]
    seur$Kit <- as.numeric(as.factor(medical_metadata$fastq[match(seur$donor,medical_metadata$donor)]))
    seur$Kit[seur$Kit!=1] <- 2
    seur$Kit[is.na(seur$Kit)] <- 1
    seur$Kit <- as.factor(seur$Kit)


    seur$mRS_3class <- as.vector(seur$mRS)
    seur$mRS_3class[seur$mRS%in% c(1,2)] <- "Good"
    seur$mRS_3class[seur$mRS%in% c(3,4)] <- "Mid"
    seur$mRS_3class[seur$mRS%in% c(5,6)] <- "Bad"
    seur$mRS_3class <- ordered(seur$mRS_3class,levels=c("Good","Mid","Bad"))
    seur$mRS_binned <- as.vector(seur$mRS)
    seur$mRS_binned[seur$mRS%in% c(1,2,3)] <- "G"
    seur$mRS_binned[seur$mRS%in% c(4,5,6)] <- "B"
    seur$mRS_binned <- ordered(seur$mRS_binned,levels=c("G","B"))
    seur$sex <- as.factor(seur$sex)
    seur$race <- as.factor(seur$race)

    saveRDS(seur,rdsoutout)
}
