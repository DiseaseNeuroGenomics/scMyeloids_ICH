# ============================================================================
# 00_setup.R
# Project: Decoding Immune Landscapes and Therapeutic Targets in Intracerebral Hemorrhage
# Author:  Dimitrios Kyriakis
# ============================================================================

# ==============================================================================
# Shared setup — source this first in every Figure1_*.R / Figure2_*.R script.
# Paths, libraries, seed, and generic save/read helpers used across all scripts.
# ==============================================================================

BASE_DIR     <- "/sc/arion/projects/CommonMind/kyriad02/Lab_Projects/ICH_Stroke/2026_08_Figure_Codes"
FIG_DIR      <- file.path(BASE_DIR, "2026_07_Paper_Figures")
TABLE_DIR    <- file.path(BASE_DIR, "2026_07_Paper_Tables")
DERIVED_DIR  <- file.path(BASE_DIR, "2026_07_Derived_Data")   # intermediate .rds objects passed between scripts
GMT_DIR      <- "/sc/arion/projects/CommonMind/kyriad02/GMT_Files/msigdb_v7.5.1_files_to_download_locally/msigdb_v7.5.1_GMTs"

dir.create(FIG_DIR,     showWarnings = FALSE, recursive = TRUE)
dir.create(TABLE_DIR,   showWarnings = FALSE, recursive = TRUE)
dir.create(DERIVED_DIR, showWarnings = FALSE, recursive = TRUE)

# ------------------------------------------------------------------------------
# Resolve this script's own directory and make it the working directory, so
# that sibling source("00_helpers_stats.R")/source("FigureN_*.R") calls resolve
# no matter where `Rscript` was invoked from. Data I/O elsewhere always uses
# file.path(BASE_DIR, ...), which is unaffected by this. Falls back to leaving
# the cwd untouched if not run via `Rscript` (e.g. sourced interactively).
# ------------------------------------------------------------------------------
.args <- commandArgs(trailingOnly = FALSE)
.file_arg <- sub("^--file=", "", .args[grepl("^--file=", .args)])
if (length(.file_arg) == 1) {
  setwd(dirname(normalizePath(.file_arg)))
}
rm(.args, .file_arg)

options(future.globals.maxSize = 100000000 * 1024^2)
set.seed(12345)

suppressPackageStartupMessages({
  # Single-cell
  library(Seurat)
  library(SingleCellExperiment)
  library(archive)
  library(scCustomize)
  library(muscat)
  library(dreamlet)
  # Visualization
  library(ComplexHeatmap)
  library(circlize)
  library(RColorBrewer)
  library(ggplot2)
  library(ggtree)
  library(cowplot)
  library(patchwork)
  library(showtext)
  # Data manipulation
  library(dplyr)
  library(tidyr)
  library(readr)
  library(broom)
  # Statistics / enrichment
  library(variancePartition)
  library(metafor)
  library(crumblr)
  library(clusterProfiler)
})

# ------------------------------------------------------------------------------
# Compressed RDS read/write — used for the large Metacells object (Figure 2)
# ------------------------------------------------------------------------------
saveCRDS <- function(object, filename, filter = "zstd") {
  stopifnot(filter %in% c("zstd", "lz4"))
  con <- archive::file_write(file = filename, filter = filter)
  open(con)
  saveRDS(object, con)
  close(con)
}

readCRDS <- function(filename) {
  con <- archive::file_read(file = filename)
  res <- readRDS(con)
  close(con)
  res
}

# ------------------------------------------------------------------------------
# Save-if-exists helpers — every panel script ends by calling these
# ------------------------------------------------------------------------------
save_plot_safe <- function(plot, filename, width, height, dpi = 600) {
  showtext_opts(dpi = dpi)
  ggsave(file.path(FIG_DIR, paste0(filename, ".pdf")), plot,
         width = width, height = height, units = "in", device = cairo_pdf, dpi = dpi)
  ggsave(file.path(FIG_DIR, paste0(filename, ".png")), plot,
         width = width, height = height, units = "in", dpi = dpi)
  cat("Saved plot:", filename, "(pdf + png)\n")
}

save_table_safe <- function(obj, filename) {
  if (!is.data.frame(obj)) obj <- as.data.frame(obj)
  write.csv(obj, file.path(TABLE_DIR, filename), row.names = FALSE)
  cat("Saved table:", filename, "(", nrow(obj), "rows,", ncol(obj), "cols )\n")
}
