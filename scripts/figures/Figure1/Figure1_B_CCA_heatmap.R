# ============================================================================
# Figure1_B_CCA_heatmap.R
# Project: Decoding Immune Landscapes and Therapeutic Targets in Intracerebral Hemorrhage
# Author:  Dimitrios Kyriakis
# ============================================================================

# ==============================================================================
# Figure 1, Panel B — Canonical Correlation Analysis (CCA) of patient metadata
#
# Pairwise canonical correlations across all clinical/demographic variables,
# to identify which metadata associate with each other and with the mRS
# outcome measures.
#
# Input:  BASE_DIR/2024_11_metadata_new.tsv
# Output: TABLE_DIR/Figure1_B_CCA_metadata.csv
#         FIG_DIR/Figure1_B_CCA_heatmap.pdf
# ==============================================================================

source("../../setup/00_setup.R")
library(pheatmap)

metadata_new <- readr::read_tsv(file.path(BASE_DIR, "2024_11_metadata_new.tsv"), col_names = TRUE)

formula <- paste0(
  " ~ mRS_binned + mRS + mRS_3class + TSH + age + dx + PHV + PEV + sex + ",
  "race + BL + CL + ICH_Score + DfCtB + NIHSS + ",
  "Suspected.Etiology..CAA + Suspected.Etiology..HTN "
)

var_part_cor <- variancePartition::canCorPairs(formula, data = metadata_new)

colnames(var_part_cor)[colnames(var_part_cor) == "Suspected.Etiology..CAA"] <- "Susp.Etiol.CAA"
colnames(var_part_cor)[colnames(var_part_cor) == "Suspected.Etiology..HTN"] <- "Susp.Etiol.HTN"
rownames(var_part_cor)[rownames(var_part_cor) == "Suspected.Etiology..CAA"] <- "Susp.Etiol.CAA"
rownames(var_part_cor)[rownames(var_part_cor) == "Suspected.Etiology..HTN"] <- "Susp.Etiol.HTN"

cca_df <- as.data.frame(var_part_cor)
cca_df$RowNames <- rownames(cca_df)

save_table_safe(cca_df, "Figure1_B_CCA_metadata.csv")

pdf(file.path(FIG_DIR, "Figure1_B_CCA_heatmap.pdf"), width = 10, height = 10)
pheatmap::pheatmap(
  var_part_cor, main = "CCA metadata", cexRow = 0.4, cexCol = 1,
  cellwidth = 25, cellheight = 25, fontsize = 15, color = brewer.pal(9, "Reds")
)
dev.off()
cat("Saved: Figure1_B_CCA_heatmap.pdf\n")
