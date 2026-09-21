#!/usr/bin/env bash
# ============================================================================
# check_pipeline_status.sh
# Project: Decoding Immune Landscapes and Therapeutic Targets in Intracerebral Hemorrhage
# Author:  Dimitrios Kyriakis
# ============================================================================

# ==============================================================================
# check_pipeline_status.sh — reports which Figure*.R scripts have already
# completed (their known final output file exists) vs still need to run, in
# recommended pipeline order, and what each one depends on.
#
# Status values:
#   DONE     - marker output file already exists
#   RUNNING  - matches a live `Rscript <name>` process (nohup) or bjobs entry
#   MISSING  - a required raw/external input file isn't present under BASE_DIR
#              (listed explicitly, e.g. "MISSING: 2024_11_Final_Immune.rds")
#   BLOCKED  - a script it depends on hasn't finished yet
#   TODO     - ready to run right now (all dependencies DONE, inputs present)
#   N/A      - console-only output, can't be checked by file existence
#
# Usage: ./check_pipeline_status.sh   (run from within R_scripts/)
# ==============================================================================

set -uo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

# Pull BASE_DIR straight out of 00_setup.R so this never drifts out of sync.
BASE_DIR=$(grep -m1 '^BASE_DIR' 00_setup.R | sed -E 's/^BASE_DIR *<- *"([^"]+)".*/\1/')
if [[ -z "$BASE_DIR" ]]; then
  echo "ERROR: could not parse BASE_DIR from 00_setup.R" >&2
  exit 1
fi
FIG_DIR="$BASE_DIR/2026_07_Paper_Figures"
TABLE_DIR="$BASE_DIR/2026_07_Paper_Tables"
DERIVED_DIR="$BASE_DIR/2026_07_Derived_Data"

# ------------------------------------------------------------------------------
# Three parallel indexed arrays (kept portable — no associative arrays, so
# this runs the same under old bash too): script name, its final-save marker
# file (or __CONSOLE_ONLY__), and a space-separated list of script names it
# depends on (pulled from each script's own "# Input:" header comment).
# ------------------------------------------------------------------------------
SCRIPTS=(
  "Figure1_00_build_pseudobulk.R"
  "Figure1_B_CCA_heatmap.R"
  "Figure1_C_umap.R"
  "Figure1_D_mRS_effect_forest.R"
  "Figure1_master_DEGs.R"
  "Figure1_E_deg_summary_bar.R"
  "Figure1_F_myeloid_enrichment.R"
  "Figure1_supp_dotplot_consistent_genes.R"
  "Figure1_analysis_QC_crosscheck.R"
  "Figure1_combined_panel.R"
  "Figure2_A_umap_metacells.R"
  "Figure2_master_DEGs.R"
  "Figure2_C_mRS_effect_tree.R"
  "Figure2_D_heatmap_top_genes.R"
  "Figure2_E_compareCluster_enrichment.R"
  "Figure2_F_00_build_SAMC_reference.R"
  "Figure2_F_SAMC_composition_bar.R"
  "Figure2_G_00_build_ExtendedGenes_loom.R"
  "Figure2_G_00_build_DEG_TFs_ExtGenes.R"
  "Figure2_G_TF_regulon_dotplot_ExtGenes.R"
  "Figure2_G_02_network_degree.R"
  "Figure2_H_hallmark_gsea_mTORC1.R"
  "Figure2_supp_enrichment_gsea_gobp_gomf.R"
  "Figure2_analysis_QC_crosscheck.R"
  "Figure3_A_irea_export.R"
  "Figure3_A_irea_cytokine_plot.R"
  "SFigure1_donor_annotations.R"
)

MARKERS=(
  "$DERIVED_DIR/Figure1_res_proc.rds"
  "$TABLE_DIR/Figure1_B_CCA_metadata.csv"
  "$FIG_DIR/Figure1_C_umap.png"
  "$FIG_DIR/Figure1_D_mRS_effect_forest.png"
  "$TABLE_DIR/Figure1_master_DEGs_all_assays.csv"
  "$FIG_DIR/Figure1_E_deg_summary_bar.png"
  "$FIG_DIR/Figure1_F_myeloid_ora_gobp.png"
  "$FIG_DIR/Figure1_supp_dotplot_consistent_genes.png"
  "$DERIVED_DIR/Figure1_leading_edge_list.rds"
  "$FIG_DIR/Figure1_combined_CDEF.png"
  "$DERIVED_DIR/Figure2_Metacells.rds"
  "$TABLE_DIR/Figure2_master_clean_DEGs_MTC456_vs_MTC123.csv"
  "$FIG_DIR/Figure2_C_composite_tree_coef.png"
  "$FIG_DIR/Figure2_D_heatmap.png"
  "$FIG_DIR/Figure2_E_GO_BP_dotplot.png"
  "$DERIVED_DIR/Figure2_F_SAMC_reference.rds"
  "$TABLE_DIR/Figure2_F_predicted_celltype_composition.csv"
  "$DERIVED_DIR/2026_08_ExtGenes_input.loom"
  "$DERIVED_DIR/Figure2_G_DEG_TFs_ExtGenes.csv"
  "$FIG_DIR/Figure2_G_TF_regulon_dotplot_ExtGenes.png"
  "$TABLE_DIR/Figure2_G_network_degree_ranking.csv"
  "$FIG_DIR/Figure2_H_hallmark_gsea_mTORC1.png"
  "$FIG_DIR/Figure2_supp_GSEA_GOMF.png"
  "__CONSOLE_ONLY__"
  "$TABLE_DIR/Figure3_IREA_ready_MTC456_vs_MTC123.xlsx"
  "$FIG_DIR/Figure3_A_irea_cytokine_plot.png"
  "$FIG_DIR/SFigure1_B_demographic_pies.png"
)

# Empty string = root script (only needs raw/external data, no prior script).
# Must stay index-aligned with SCRIPTS/MARKERS above.
DEPENDS=(
  ""                                                                                  # Figure1_00_build_pseudobulk.R
  ""                                                                                  # Figure1_B_CCA_heatmap.R
  ""                                                                                  # Figure1_C_umap.R
  "Figure1_00_build_pseudobulk.R"                                                    # Figure1_D_mRS_effect_forest.R
  "Figure1_00_build_pseudobulk.R"                                                    # Figure1_master_DEGs.R
  "Figure1_master_DEGs.R"                                                            # Figure1_E_deg_summary_bar.R
  "Figure1_master_DEGs.R"                                                            # Figure1_F_myeloid_enrichment.R
  "Figure1_master_DEGs.R"                                                            # Figure1_supp_dotplot_consistent_genes.R
  "Figure1_master_DEGs.R Figure1_F_myeloid_enrichment.R"                             # Figure1_analysis_QC_crosscheck.R
  "Figure1_C_umap.R Figure1_D_mRS_effect_forest.R Figure1_E_deg_summary_bar.R Figure1_F_myeloid_enrichment.R"  # Figure1_combined_panel.R
  ""                                                                                  # Figure2_A_umap_metacells.R
  "Figure2_A_umap_metacells.R"                                                       # Figure2_master_DEGs.R
  "Figure2_A_umap_metacells.R"                                                       # Figure2_C_mRS_effect_tree.R
  "Figure2_A_umap_metacells.R Figure2_master_DEGs.R"                                 # Figure2_D_heatmap_top_genes.R
  "Figure2_A_umap_metacells.R"                                                       # Figure2_E_compareCluster_enrichment.R
  ""                                                                                  # Figure2_F_00_build_SAMC_reference.R
  "Figure2_F_00_build_SAMC_reference.R Figure2_A_umap_metacells.R"                   # Figure2_F_SAMC_composition_bar.R
  "Figure2_A_umap_metacells.R"                                                       # Figure2_G_00_build_ExtendedGenes_loom.R
  "Figure2_A_umap_metacells.R"                                                       # Figure2_G_00_build_DEG_TFs_ExtGenes.R
  "Figure2_G_00_build_DEG_TFs_ExtGenes.R Figure2_master_DEGs.R"                      # Figure2_G_TF_regulon_dotplot_ExtGenes.R
  "Figure2_G_00_build_DEG_TFs_ExtGenes.R"                                            # Figure2_G_02_network_degree.R
  "Figure2_master_DEGs.R"                                                            # Figure2_H_hallmark_gsea_mTORC1.R
  "Figure2_master_DEGs.R"                                                            # Figure2_supp_enrichment_gsea_gobp_gomf.R
  "Figure1_analysis_QC_crosscheck.R Figure2_master_DEGs.R Figure1_master_DEGs.R"     # Figure2_analysis_QC_crosscheck.R
  "Figure2_master_DEGs.R"                                                            # Figure3_A_irea_export.R
  "Figure3_A_irea_export.R"                                                          # Figure3_A_irea_cytokine_plot.R
  ""                                                                                  # SFigure1_donor_annotations.R
)

# Raw/external files this script needs that no other script in this pipeline
# produces (root scripts' source data, plus a couple of externally-supplied
# files). Space-separated full paths; empty = nothing beyond what DEPENDS
# already covers. Every entry here should live under BASE_DIR (per current
# convention: everything consolidated into the flat working folder) — full
# paths so this doubles as "where it's expected."
RAW_INPUTS=(
  "$BASE_DIR/2024_11_Final_Immune.rds"                                                                                # Figure1_00_build_pseudobulk.R
  "$BASE_DIR/2024_11_metadata_new.tsv"                                                                                # Figure1_B_CCA_heatmap.R
  "$BASE_DIR/2024_11_Final_Immune.rds"                                                                                # Figure1_C_umap.R
  ""                                                                                                                   # Figure1_D_mRS_effect_forest.R
  ""                                                                                                                   # Figure1_master_DEGs.R
  ""                                                                                                                   # Figure1_E_deg_summary_bar.R
  ""                                                                                                                   # Figure1_F_myeloid_enrichment.R
  ""                                                                                                                   # Figure1_supp_dotplot_consistent_genes.R
  ""                                                                                                                   # Figure1_analysis_QC_crosscheck.R
  ""                                                                                                                   # Figure1_combined_panel.R
  "$BASE_DIR/2024_03_28_Myeloid_Metacells_Subclass_ADAM.rds.ztsd"                                                     # Figure2_A_umap_metacells.R
  ""                                                                                                                   # Figure2_master_DEGs.R
  ""                                                                                                                   # Figure2_C_mRS_effect_tree.R
  ""                                                                                                                   # Figure2_D_heatmap_top_genes.R
  ""                                                                                                                   # Figure2_E_compareCluster_enrichment.R
  "$BASE_DIR/External_Data/GSE189432_annotations.csv.gz"                                                              # Figure2_F_00_build_SAMC_reference.R (canary file only — doesn't check every GSM sample triplet)
  ""                                                                                                                   # Figure2_F_SAMC_composition_bar.R
  "$BASE_DIR/database/allTFs_hg38.txt $BASE_DIR/database/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl"                  # Figure2_G_00_build_ExtendedGenes_loom.R
  "$DERIVED_DIR/2026_08_ExtGenes_RegulonsAUC.csv"                                                                     # Figure2_G_00_build_DEG_TFs_ExtGenes.R (pipelines/extract_loom_regulons.py, run manually after the bsub GRN job)
  "$DERIVED_DIR/2026_08_ExtGenes_Regulons.csv"                                                                        # Figure2_G_TF_regulon_dotplot_ExtGenes.R (same extract_loom_regulons.py run)
  "$DERIVED_DIR/2026_08_ExtGenes_Regulons.csv $DERIVED_DIR/2026_08_ExtGenes_adj.csv"                                  # Figure2_G_02_network_degree.R
  ""                                                                                                                   # Figure2_H_hallmark_gsea_mTORC1.R
  ""                                                                                                                   # Figure2_supp_enrichment_gsea_gobp_gomf.R
  ""                                                                                                                   # Figure2_analysis_QC_crosscheck.R
  ""                                                                                                                   # Figure3_A_irea_export.R
  "$TABLE_DIR/Figure3_IREA_Results.csv"                                                                               # Figure3_A_irea_cytokine_plot.R (manually downloaded from IREA)
  "$BASE_DIR/2024_11_Final_Immune.rds"                                                                                # SFigure1_donor_annotations.R
)

n=${#SCRIPTS[@]}

# ------------------------------------------------------------------------------
# Is a given script currently running, as a local nohup job or an LSF job?
# ------------------------------------------------------------------------------
is_running() {
  local script="$1"
  if pgrep -f "Rscript[ ]+${script}\b" >/dev/null 2>&1; then
    return 0
  fi
  if command -v bjobs >/dev/null 2>&1; then
    if bjobs -w 2>/dev/null | grep -q "$script"; then
      return 0
    fi
  fi
  return 1
}

index_of() {
  local target="$1" i
  for ((i = 0; i < n; i++)); do
    [[ "${SCRIPTS[$i]}" == "$target" ]] && { echo "$i"; return 0; }
  done
  echo "-1"
}

# ------------------------------------------------------------------------------
# Pass 1: base status per script (DONE / RUNNING / TODO / N/A), ignoring
# dependencies and raw inputs for now.
# ------------------------------------------------------------------------------
STATUSES=()
for ((i = 0; i < n; i++)); do
  marker="${MARKERS[$i]}"
  script="${SCRIPTS[$i]}"
  if is_running "$script"; then
    STATUSES+=("RUNNING")
  elif [[ "$marker" == "__CONSOLE_ONLY__" ]]; then
    STATUSES+=("N/A")
  elif [[ -f "$marker" ]]; then
    STATUSES+=("DONE")
  else
    STATUSES+=("TODO")
  fi
done

# ------------------------------------------------------------------------------
# Pass 2: demote TODO -> MISSING if a required raw/external input file isn't
# present under BASE_DIR. Records which specific file(s) for display.
# ------------------------------------------------------------------------------
MISSING_FILES=()
for ((i = 0; i < n; i++)); do
  MISSING_FILES+=("")
  [[ "${STATUSES[$i]}" != "TODO" ]] && continue
  missing=""
  for f in ${RAW_INPUTS[$i]}; do
    [[ -f "$f" ]] || missing="$missing $(basename "$f")"
  done
  if [[ -n "$missing" ]]; then
    STATUSES[$i]="MISSING"
    MISSING_FILES[$i]="${missing# }"
  fi
done

# ------------------------------------------------------------------------------
# Pass 3: demote remaining TODO -> BLOCKED if any dependency isn't DONE yet.
# ------------------------------------------------------------------------------
for ((i = 0; i < n; i++)); do
  [[ "${STATUSES[$i]}" != "TODO" ]] && continue
  for dep in ${DEPENDS[$i]}; do
    dep_idx=$(index_of "$dep")
    [[ "$dep_idx" == "-1" ]] && continue
    if [[ "${STATUSES[$dep_idx]}" != "DONE" ]]; then
      STATUSES[$i]="BLOCKED"
      break
    fi
  done
done

# ------------------------------------------------------------------------------
# Color helpers (skip if not a terminal)
# ------------------------------------------------------------------------------
if [[ -t 1 ]]; then
  GREEN=$'\033[32m'; YELLOW=$'\033[33m'; RED=$'\033[31m'; ORANGE=$'\033[33m'; DIM=$'\033[2m'; RESET=$'\033[0m'
else
  GREEN=""; YELLOW=""; RED=""; ORANGE=""; DIM=""; RESET=""
fi

printf "%-42s %-9s %-55s %s\n" "SCRIPT" "STATUS" "DEPENDS ON" "MARKER"
printf '%s\n' "-------------------------------------------------------------------------------------------------------------------"

n_done=0; n_running=0; n_todo=0; n_blocked=0; n_na=0; n_missing=0
for ((i = 0; i < n; i++)); do
  script="${SCRIPTS[$i]}"
  marker="${MARKERS[$i]}"
  status="${STATUSES[$i]}"
  deps="${DEPENDS[$i]}"
  deps_display="${deps// /, }"
  [[ -z "$deps_display" ]] && deps_display="-"

  case "$status" in
    DONE)    color="$GREEN";  n_done=$((n_done+1)) ;;
    RUNNING) color="$YELLOW"; n_running=$((n_running+1)) ;;
    MISSING) color="$RED";    n_missing=$((n_missing+1)) ;;
    BLOCKED) color="$RED";    n_blocked=$((n_blocked+1)) ;;
    N/A)     color="$DIM";    n_na=$((n_na+1)) ;;
    *)       color="";        n_todo=$((n_todo+1)) ;;
  esac

  if [[ "$status" == "MISSING" ]]; then
    marker_display="MISSING FILE: ${MISSING_FILES[$i]}"
  elif [[ "$marker" == "__CONSOLE_ONLY__" ]]; then
    marker_display="(console-only - check its log file)"
  else
    marker_display="$(basename "$marker")"
  fi

  printf "%-42s ${color}%-9s${RESET} %-55s %s\n" "$script" "$status" "$deps_display" "$marker_display"
done

echo
echo "Done: $n_done   Running: $n_running   Ready (TODO): $n_todo   Blocked: $n_blocked   Missing input: $n_missing   N/A: $n_na"
