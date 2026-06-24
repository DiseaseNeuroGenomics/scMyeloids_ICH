#!/usr/bin/env python3
# Author: Dimitrios Kyriakis
# Description: scDRS polygenic enrichment scoring at single-cell level.
#              Loads myeloid metacell data, scores cells for each trait in
#              the provided gene-set file, and saves per-trait score and
#              group-level summary CSVs.
# Usage: python 8.scDRS.py
#   --h5ad         Path to myeloid metacell h5ad
#   --gs_file      Path to gene-set file (.gs format)
#   --output_dir   Directory for output CSVs (default: results/2024_04_01_scDRS/)
#   --group_col    Metadata column for group analysis (default: Metacells_identity)

import argparse
import os
import warnings

import numpy as np
import pandas as pd
import scanpy as sc
import scdrs

warnings.filterwarnings("ignore")

parser = argparse.ArgumentParser(description="scDRS polygenic enrichment scoring")
parser.add_argument("--h5ad",       required=True,
                    help="Path to myeloid metacell h5ad")
parser.add_argument("--gs_file",    required=True,
                    help="Gene-set file (.gs) with trait definitions")
parser.add_argument("--output_dir", default="results/2024_04_01_scDRS",
                    help="Output directory for CSVs")
parser.add_argument("--group_col",  default="Metacells_identity",
                    help="Metadata column for group-level analysis")
args = parser.parse_args()

import random
random.seed(12345678)
np.random.seed(12345678)

os.makedirs(args.output_dir, exist_ok=True)
traits_dir  = os.path.join(args.output_dir, "traits")
summary_dir = os.path.join(args.output_dir, "traits", "summary")
os.makedirs(traits_dir, exist_ok=True)
os.makedirs(summary_dir, exist_ok=True)

# ── Load data ─────────────────────────────────────────────────────────────────
adata = sc.read_h5ad(args.h5ad)
adata.X = adata.layers.get("data", adata.X)

# ── Load gene sets ────────────────────────────────────────────────────────────
dict_gs = scdrs.util.load_gs(
    args.gs_file,
    src_species="human",
    dst_species="human",
    to_intersect=adata.var_names
)

# ── Define trait groups ───────────────────────────────────────────────────────
all_traits = list(dict_gs.keys())
print(f"Loaded {len(all_traits)} traits from {args.gs_file}")

# ── Pre-process ───────────────────────────────────────────────────────────────
scdrs.preprocess(adata, n_mean_bin=20, n_var_bin=20, copy=False)
sc.pp.neighbors(adata)

# ── Score each trait ──────────────────────────────────────────────────────────
for trait in all_traits:
    summary_path = os.path.join(summary_dir,
                                f"Metacells_identity_{trait}_summary.csv")
    score_path   = os.path.join(traits_dir,
                                f"Metacells_identity_{trait}.csv")

    if os.path.exists(summary_path):
        print(f"Skipping {trait} — summary exists")
        continue

    gene_list, gene_weights = dict_gs[trait]
    df_score = scdrs.score_cell(
        data=adata,
        gene_list=gene_list,
        gene_weight=gene_weights,
        ctrl_match_key="mean_var",
        n_ctrl=100,
        weight_opt="vs",
        return_ctrl_raw_score=False,
        return_ctrl_norm_score=True,
        verbose=False
    )

    df_score.to_csv(score_path, sep="\t", escapechar=" ")

    df_stats = scdrs.method.downstream_group_analysis(
        adata=adata,
        df_full_score=df_score,
        group_cols=[args.group_col]
    )[args.group_col]
    df_stats.to_csv(summary_path, sep="\t", escapechar=" ")

    print(f"Saved trait: {trait}")

print(f"Done. Results in: {args.output_dir}")
