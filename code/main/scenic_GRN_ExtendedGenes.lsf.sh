#!/bin/bash
# ============================================================================
# scenic_GRN_ExtendedGenes.lsf.sh
# Project: Decoding Immune Landscapes and Therapeutic Targets in Intracerebral Hemorrhage
# Author:  Dimitrios Kyriakis
# ============================================================================

# ==============================================================================
# pySCENIC GRN pipeline, run against 2026_08_ExtGenes_input.loom — a rebuild
# of the original ExtendedGenes gene-curation logic
# (4.Figure2_G_Prepare_For_SCENIC.ipynb), produced by
# Rscripts/Figure2_G_00_build_ExtendedGenes_loom.R. Run that R script first.
#
# Fully self-contained: every input this job reads lives inside
# 2026_08_Figure_Codes (database/, pipelines/, 2026_07_Derived_Data), and all
# writes stay confined to 2026_07_Derived_Data. Nothing outside this project
# folder is read or written.
#
# Adapted from pipelines/scenic_GRN_ICH.lsf.sh (the run that produced
# 2025_MetaCL_auc_v10.loom from a DIFFERENT input loom,
# 2024_04_01_Myeloid_Metacells_pcaENTREZids.loom, kept as-is and untouched —
# this is a genuinely separate GRN inference run, new output names throughout
# so the two can be compared side by side rather than one overwriting the other).
#
# Only computes the v10 motif database (Methods section confirms pySCENIC was
# run with v10).
#
# Runs 3 steps, each skipped if its output already exists:
#   1. GRNBoost2 (adjacencies)     -> OUTPUT_ADJ
#   2. pyscenic ctx, motif db v10  -> OUTPUT_REG_v10
#   3. pyscenic aucell             -> OUTPUT_AUCELL_v10
#
# Prerequisite: Rscript Rscripts/Figure2_G_00_build_ExtendedGenes_loom.R
#               (produces 2026_07_Derived_Data/2026_08_ExtGenes_input.loom)
# ==============================================================================

#BSUB -J GRN_ExtendedGenes
#BSUB -P acc_CommonMind
#BSUB -q premium
#BSUB -n 40
#BSUB -M 80G
#BSUB -W 90:00
#BSUB -R span[hosts=1]
#BSUB -o /sc/arion/projects/CommonMind/kyriad02/Lab_Projects/ICH_Stroke/2026_08_Figure_Codes/2026_07_Derived_Data/2026_08_GRN_ExtendedGenes.out
#BSUB -oo /sc/arion/projects/CommonMind/kyriad02/Lab_Projects/ICH_Stroke/2026_08_Figure_Codes/2026_07_Derived_Data/2026_08_GRN_ExtendedGenes.stdout
#BSUB -eo /sc/arion/projects/CommonMind/kyriad02/Lab_Projects/ICH_Stroke/2026_08_Figure_Codes/2026_07_Derived_Data/2026_08_GRN_ExtendedGenes.stderr
#BMAIL -u  kyriakds@gmail.com
#BMAIL -N
#BMAIL -B

ml purge
ml gcc/8.3.0
ml cmake/3.22.0
ml libpng
ml libxml2
ml libxmlsec1
ml gdal
ml hdf5/1.12.1
ml singularity

FIGURE_CODES_DIR="/sc/arion/projects/CommonMind/kyriad02/Lab_Projects/ICH_Stroke/2026_08_Figure_Codes"
WORK_OUT_DIR="$FIGURE_CODES_DIR/2026_07_Derived_Data"
mkdir -p "$WORK_OUT_DIR"
cd "$WORK_OUT_DIR"

# ---- Inputs — all inside 2026_08_Figure_Codes ----
LOOM_FILE="$WORK_OUT_DIR/2026_08_ExtGenes_input.loom"
ALBERTO_SCRIPT="$FIGURE_CODES_DIR/pipelines/alberto_multi.py"
TFS_FILE="$FIGURE_CODES_DIR/database/allTFs_hg38.txt"
RANKINGS_FILE_v10="$FIGURE_CODES_DIR/database/hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather"
ANNOTATIONS_FILE_v10="$FIGURE_CODES_DIR/database/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl"
SIF_IMAGE="$FIGURE_CODES_DIR/database/aertslab-pyscenic-0.12.1.sif"

if [ ! -e "$LOOM_FILE" ]; then
    echo "ERROR: $LOOM_FILE not found. Run Rscripts/Figure2_G_00_build_ExtendedGenes_loom.R first."
    exit 1
fi

# ---- Writes — all confined to WORK_OUT_DIR ----
OUTPUT_ADJ="$WORK_OUT_DIR/2026_08_ExtGenes_adj.csv"
OUTPUT_REG_v10="$WORK_OUT_DIR/2026_08_ExtGenes_reg_v10.csv"
OUTPUT_AUCELL_v10="$WORK_OUT_DIR/2026_08_ExtGenes_auc_v10.loom"

METHOD="grnboost2"
NUM_WORKERS=40
SEED=777
MASK_DROPOUTS="--mask_dropouts"

# Only need to bind FIGURE_CODES_DIR now — every input/output lives under it.
SINGULARITY_BIND="--bind $FIGURE_CODES_DIR:$FIGURE_CODES_DIR"

if [ ! -e "$OUTPUT_ADJ" ]; then
    singularity run $SINGULARITY_BIND \
        "$SIF_IMAGE" \
        python "$ALBERTO_SCRIPT" \
        "$LOOM_FILE" \
        "$TFS_FILE" \
        --method "$METHOD" \
        --output "$OUTPUT_ADJ" \
        --num_workers "$NUM_WORKERS" \
        --seed "$SEED"
else
    echo "Skipping execution: '$OUTPUT_ADJ' already exists."
fi

if [ ! -e "$OUTPUT_REG_v10" ]; then
    singularity run $SINGULARITY_BIND \
        "$SIF_IMAGE" \
        pyscenic ctx \
            "$OUTPUT_ADJ" \
            "$RANKINGS_FILE_v10" \
            --annotations_fname "$ANNOTATIONS_FILE_v10" \
            --expression_mtx_fname "$LOOM_FILE" \
            --output "$OUTPUT_REG_v10" \
            $MASK_DROPOUTS \
            --num_workers "$NUM_WORKERS"
else
    echo "Skipping execution: '$OUTPUT_REG_v10' already exists."
fi

if [ ! -e "$OUTPUT_AUCELL_v10" ]; then
    singularity run $SINGULARITY_BIND \
        "$SIF_IMAGE" \
        pyscenic aucell \
            "$LOOM_FILE" \
            "$OUTPUT_REG_v10" \
            --output "$OUTPUT_AUCELL_v10" \
            --num_workers "$NUM_WORKERS"
else
    echo "Skipping execution: '$OUTPUT_AUCELL_v10' already exists."
fi
