#!/usr/bin/env Rscript
# Author: Dimitrios Kyriakis
# Description: Gene Ontology (Molecular Function) enrichment analysis between
#              MTC_123 and MTC_456 metacluster DEGs using clusterProfiler.
#              Produces SFigure 3B dotplot.
# Usage: Rscript SFigure_3B.R <metacell_rds> <metadata_csv> <output_pdf>
#   <metacell_rds>   Path to 2024_03_28_Myeloid_Metacells_Subclass_ADAM.rds.ztsd
#   <metadata_csv>   Path to Data/2024_Metadata_ICH.csv
#   <output_pdf>     Output PDF path (default: SFigure_3B.pdf)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
    stop("Usage: Rscript SFigure_3B.R <metacell_rds> <metadata_csv> [output_pdf]")
}

input_rds    <- args[1]
metadata_csv <- args[2]
output_pdf   <- if (length(args) >= 3) args[3] else "SFigure_3B.pdf"

set.seed(12345)

suppressPackageStartupMessages({
    library(Seurat)
    library(clusterProfiler)
    library(org.Hs.eg.db)
    library(ggplot2)
    library(dplyr)
    library(readr)
})

source('utils.R')

Metacells <- readCRDS(input_rds)

New_Metadata <- read_csv(metadata_csv, show_col_types = FALSE)
rownames(New_Metadata) <- New_Metadata$donor

Metacells$TwoClusters <- as.vector(Metacells$Clusters)
Metacells$TwoClusters[Metacells$Clusters %in% c("MTC_1", "MTC_2", "MTC_3")] <- "MTC_123"
Metacells$TwoClusters[Metacells$Clusters %in% c("MTC_4", "MTC_5", "MTC_6")] <- "MTC_456"

if ("ICH Score" %in% colnames(New_Metadata)) {
    Metacells$ICH         <- as.numeric(unlist(New_Metadata[Metacells$Donor, "ICH Score"]))
    Metacells$mRS_binned  <- ifelse(Metacells$mRS >= 4, "B", "G")
}

Metacells_2 <- subset(Metacells, subset = TwoClusters %in% c("MTC_123", "MTC_456"))
Idents(Metacells_2) <- "TwoClusters"

deg_pa <- FindAllMarkers(Metacells_2, only.pos = TRUE)
significant_genes <- deg_pa %>% filter(p_val_adj < 0.01, avg_log2FC > 1.5)

formula_MF <- compareCluster(
    gene ~ cluster,
    data           = significant_genes,
    fun            = "enrichGO",
    OrgDb          = org.Hs.eg.db,
    keyType        = "SYMBOL",
    ont            = "MF",
    pAdjustMethod  = "BH",
    pvalueCutoff   = 0.05
)

p <- dotplot(formula_MF, showCategory = 10)

pdf(output_pdf, width = 7, height = 10)
print(p)
dev.off()

message("Saved: ", output_pdf)
