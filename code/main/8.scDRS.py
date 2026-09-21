#!/usr/bin/env python3
# Author: Dimitrios Kyriakis
# Description: scDRS polygenic enrichment scoring at single-cell level.
#              Loads myeloid metacell data, scores cells for each trait in
#              the provided gene-set file, and saves per-trait score and
#              group-level summary CSVs.
# Usage: python 8.scDRS.py
#   --h5ad         Path to myeloid metacell h5ad
#   --gs_file      Path to gene-set file (.gs format)
#   --n_ctrl       Number of control gene sets per trait (default: 100).
#                  Determines the smallest attainable Monte-Carlo p-value,
#                  which is 1/(n_ctrl+1).
#   --output_dir   Directory for output CSVs. REQUIRED — there is no default,
#                  specifically so a rerun cannot silently land on top of an
#                  existing results directory. Use a fresh, empty path for
#                  every run (e.g. results/scDRS_nctrl100_recheck,
#                  results/scDRS_nctrl1000). Never point this at an existing
#                  results directory you want to keep.
#   --group_col    Metadata column for group analysis (default: Metacells_identity)
#   --resume       Opt-in: skip traits whose summary CSV already exists in
#                  output_dir, instead of erroring. Off by default. This
#                  used to be the unconditional behaviour, which is exactly
#                  what caused a later n_ctrl=1000 rerun to silently reuse
#                  n_ctrl=100 results — every trait already had a summary
#                  file on disk, so nothing was recomputed and the run
#                  finished instantly without warning.

import argparse
import os
import sys
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
parser.add_argument("--n_ctrl",     type=int, default=100,
                    help="Number of control gene sets per trait (default: 100)")
parser.add_argument("--output_dir", required=True,
                    help="Output directory for CSVs. Must not already contain "
                         "a traits/summary/ directory unless --resume is given.")
parser.add_argument("--group_col",  default="Metacells_identity",
                    help="Metadata column for group-level analysis")
parser.add_argument("--resume", action="store_true",
                    help="Skip traits whose summary CSV already exists in "
                         "output_dir instead of erroring. Use only when "
                         "deliberately continuing an interrupted run at the "
                         "SAME n_ctrl — never after changing n_ctrl.")
args = parser.parse_args()

import random
random.seed(12345678)
np.random.seed(12345678)

traits_dir  = os.path.join(args.output_dir, "traits")
summary_dir = os.path.join(args.output_dir, "traits", "summary")

existing_summaries = (
    [f for f in os.listdir(summary_dir) if f.endswith(".csv")]
    if os.path.isdir(summary_dir) else []
)
if existing_summaries and not args.resume:
    sys.exit(
        f"ERROR: {summary_dir} already contains {len(existing_summaries)} "
        f"summary file(s). Refusing to run into it, because per-trait "
        f"results would silently be skipped rather than recomputed, which "
        f"is how the n_ctrl=100 vs n_ctrl=1000 comparison ended up "
        f"identical last time.\n"
        f"Fix: pass a fresh --output_dir (e.g. "
        f"results/scDRS_nctrl{args.n_ctrl}_$(date +%Y%m%d)), or pass "
        f"--resume if you are deliberately continuing an interrupted run "
        f"at the SAME n_ctrl={args.n_ctrl}."
    )

os.makedirs(traits_dir, exist_ok=True)
os.makedirs(summary_dir, exist_ok=True)

print(f"n_ctrl={args.n_ctrl}  ->  Monte-Carlo p-value floor = "
      f"1/{args.n_ctrl + 1} = {1 / (args.n_ctrl + 1):.6f}")
print(f"Writing to: {args.output_dir}")

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
        if not args.resume:
            # Should be unreachable given the pre-flight check above, but
            # kept as a hard stop in case files appear mid-run.
            sys.exit(f"ERROR: {summary_path} exists and --resume was not set.")
        print(f"Skipping {trait} — summary exists (--resume)")
        continue

    gene_list, gene_weights = dict_gs[trait]
    df_score = scdrs.score_cell(
        data=adata,
        gene_list=gene_list,
        gene_weight=gene_weights,
        ctrl_match_key="mean_var",
        n_ctrl=args.n_ctrl,
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

print(f"Done. n_ctrl={args.n_ctrl}. Results in: {args.output_dir}")
