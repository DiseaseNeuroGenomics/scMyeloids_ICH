#!/usr/bin/env Rscript
# ============================================================================
# build_supplementary_tables.R
# Project: Decoding Immune Landscapes and Therapeutic Targets in Intracerebral Hemorrhage
# Author:  Dimitrios Kyriakis
# ============================================================================

# ==============================================================================
# Combine the Figure 4 statistics into two multi-sheet supplementary tables
# plus one Source Data workbook.
#
# Run from MOUSE_EXPERIMENTS/
# ==============================================================================

suppressPackageStartupMessages({ library(tidyverse); library(writexl) })

BEHAV <- "../tables/mouse"                      # behavioural tables
MORPH <- "../tables/mouse"                  # morphology tables
OUT   <- "../tables/supplementary"
dir.create(OUT, showWarnings = FALSE)

rd <- function(dir, f) {
  p <- file.path(dir, f)
  if (!file.exists(p)) { message("MISSING: ", p); return(NULL) }
  read.csv(p, check.names = FALSE)
}

keep_primary <- function(x) if (!is.null(x) && "Set" %in% names(x))
  filter(x, Set == "primary") %>% select(-Set) else x

# ------------------------------------------------------------------ S.Table 11
# Behavioural analysis (Figure 4c)
st11 <- list(
  `A. Per-timepoint contrasts` = keep_primary(rd(BEHAV, "03_per_timepoint_contrasts.csv")),
  `B. Trajectory contrasts`    = keep_primary(rd(BEHAV, "04_slope_contrasts.csv")),
  `C. Baseline balance`        = rd(BEHAV, "06_baseline_balance_check.csv"),
  `D. Marginal means`          = keep_primary(rd(BEHAV, "05_estimated_marginal_means.csv")),
  `E. Model coefficients`      = keep_primary(rd(BEHAV, "01_model_coefficients.csv")),
  `F. n per group per time`    = keep_primary(rd(BEHAV, "07_n_per_group_per_timepoint.csv"))
) %>% compact()
write_xlsx(st11, file.path(OUT, "SupplementaryTable11_behaviour.xlsx"))

# ------------------------------------------------------------------ S.Table 12
# Microglial morphology analysis (Figure 4d,e)
st12 <- list(
  `A. Distance model`          = rd(MORPH, "T8_morphology_distance_model_coefficients.csv"),
  `B. Predicted probabilities` = rd(MORPH, "T9_morphology_predicted_probabilities.csv"),
  `C. Windowed model`          = rd(MORPH, "T11_morphology_window_model_coefficients.csv"),
  `D. Cohort counts`           = rd(MORPH, "T13_morphology_cohort_counts.csv"),
  `E. Cells per morphotype`    = rd(MORPH, "T14_morphology_cell_counts_by_morphotype.csv"),
  `F. Animal-level slopes`     = rd(MORPH, "T15_morphology_animal_level_slopes.csv"),
  `G. Animal-level contrasts`  = rd(MORPH, "T16_morphology_animal_level_contrasts.csv")
) %>% compact()
write_xlsx(st12, file.path(OUT, "SupplementaryTable12_morphology.xlsx"))

# ------------------------------------------------------------------ Source Data
# Nature requires the underlying values for every plotted graph.
src <- list(
  `Fig4c per-animal`   = rd(MORPH, "T7_behaviour_source_data_per_animal.csv"),
  `Fig4c summary`      = keep_primary(rd(BEHAV, "08_descriptives_source_data.csv")),
  `Fig4d probabilities`= rd(MORPH, "T9_morphology_predicted_probabilities.csv"),
  `Fig4e per-animal`   = rd(MORPH, "T12_morphology_source_data_per_animal.csv")
) %>% compact()
write_xlsx(src, file.path(OUT, "SourceData_Figure4.xlsx"))

message("\nWrote ", length(list.files(OUT)), " files to ", normalizePath(OUT))

src
