# ============================================================================
# Figure3_A_irea_export.R
# Project: Decoding Immune Landscapes and Therapeutic Targets in Intracerebral Hemorrhage
# Author:  Dimitrios Kyriakis
# ============================================================================

# ==============================================================================
# Figure 3, Panel A (input build) — export MTC_456 vs MTC_123 log2FC in IREA's
# required format (MaxAbs-scaled) for upload to the IREA web tool
# (https://immunet.mssm.edu/irea/). IREA was run with species=human,
# cell type=macrophage, Gene Diff Cutoff=0.07 (Methods). IREA's output CSV
# must be downloaded manually and saved to TABLE_DIR/Figure3_IREA_Results.csv
# before running Figure3_A_irea_cytokine_plot.R.
#
# Input:  TABLE_DIR/Figure2_master_clean_DEGs_MTC456_vs_MTC123.csv (Figure2_master_DEGs.R)
# Output: TABLE_DIR/Figure3_IREA_ready_MTC456_vs_MTC123.xlsx
# ==============================================================================

source("../../setup/00_setup.R")
library(data.table)
library(writexl)

clean_DEGs <- read.csv(file.path(TABLE_DIR, "Figure2_master_clean_DEGs_MTC456_vs_MTC123.csv"))

dt_irea <- as.data.table(clean_DEGs)[, .(Gene = gene, MTC456_vs_MTC123 = log2FC)]
dt_irea[, MTC456_vs_MTC123 := MTC456_vs_MTC123 / max(abs(MTC456_vs_MTC123), na.rm = TRUE)]

write_xlsx(dt_irea, file.path(TABLE_DIR, "Figure3_IREA_ready_MTC456_vs_MTC123.xlsx"))
cat("Saved: Figure3_IREA_ready_MTC456_vs_MTC123.xlsx — upload this to IREA, then save its\n",
    "output as", file.path(TABLE_DIR, "Figure3_IREA_Results.csv"), "\n")
head(dt_irea)
