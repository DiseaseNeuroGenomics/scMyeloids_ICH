#!/usr/bin/env python3
# Author: Dimitrios Kyriakis
# Description: LIANA cell-cell interaction scoring across immune metacell clusters.
#              Runs multiple CCI methods (CellPhoneDB, CellChat, NATMI, etc.)
#              and stores results in adata.uns["liana_res"].
#              Saves the updated h5ad and a flat CSV of results.
# Usage: python 9.Liana.py
#   --h5ad        Path to input immune metacell h5ad
#   --output_h5ad Path for output h5ad (default: <input>_liana.h5ad)
#   --output_csv  Path for flat CCI CSV (default: liana_results.csv)
#   --group_col   Cell-type metadata column (default: Clusters)
#   --n_jobs      Number of parallel jobs (default: 4)

import argparse
import os
import warnings

warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
import scanpy as sc
import liana as li
from liana.method import (singlecellsignalr, connectome, cellphonedb,
                           natmi, logfc, cellchat, geometric_mean)

parser = argparse.ArgumentParser(description="LIANA CCI scoring")
parser.add_argument("--h5ad",        required=True,
                    help="Input immune metacell h5ad")
parser.add_argument("--output_h5ad", default=None,
                    help="Output h5ad path")
parser.add_argument("--output_csv",  default="liana_results.csv",
                    help="Output CSV of liana_res")
parser.add_argument("--group_col",   default="Clusters",
                    help="Cell-type column in obs (default: Clusters)")
parser.add_argument("--n_jobs",      type=int, default=4,
                    help="Parallel jobs (default: 4)")
args = parser.parse_args()

np.random.seed(12345)

if args.output_h5ad is None:
    base, ext = os.path.splitext(args.h5ad)
    args.output_h5ad = base + "_liana" + ext

adata = sc.read_h5ad(args.h5ad)
print(f"Loaded: {adata.shape[0]} cells × {adata.shape[1]} genes")
print(f"Cell types in '{args.group_col}': {adata.obs[args.group_col].unique().tolist()}")

# Use raw counts for LIANA
if "counts" in adata.layers:
    adata.X = adata.layers["counts"].copy()

sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# Run LIANA consensus (all built-in methods)
li.mt.rank_aggregate(
    adata,
    groupby  = args.group_col,
    use_raw  = False,
    verbose  = True,
    n_jobs   = args.n_jobs
)

adata.write_h5ad(args.output_h5ad)
print(f"Saved h5ad: {args.output_h5ad}")

# Export flat CSV
liana_res = adata.uns["liana_res"].copy()
liana_res.to_csv(args.output_csv, index=False)
print(f"Saved CSV: {args.output_csv}")
