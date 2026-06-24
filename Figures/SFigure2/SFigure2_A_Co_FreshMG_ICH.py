#!/usr/bin/env python3
# Author: Dimitrios Kyriakis
# Description: Pearson correlation heatmap comparing FreshMG pseudobulk
#              expression profiles with ICH biopsy profiles.
#              Produces SFigure 2A.
# Usage: python SFigure2_A_Co_FreshMG_ICH.py
#   --fresh_h5ad   Path to FreshMG h5ad (pseudobulk aggregation)
#   --biopsies_h5ad  Path to ICH biopsies h5ad
#   --fresh_csv    Path to pre-computed FreshMG pseudobulk CSV (if already aggregated)
#   --biopsies_csv Path to pre-computed biopsies pseudobulk CSV (if already aggregated)
#   --output_pdf   Output PDF path (default: SFigures/SFigure2A.pdf)

import argparse
import os

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import scipy.stats as stats

parser = argparse.ArgumentParser(
    description="Pearson correlation heatmap: FreshMG vs ICH biopsies"
)
parser.add_argument("--fresh_csv",    required=True,
                    help="Pre-aggregated FreshMG pseudobulk CSV (genes x subtypes)")
parser.add_argument("--biopsies_csv", required=True,
                    help="Pre-aggregated ICH biopsies pseudobulk CSV (genes x subtypes)")
parser.add_argument("--output_pdf",   default="SFigures/SFigure2A.pdf",
                    help="Output PDF path")
args = parser.parse_args()

os.makedirs(os.path.dirname(args.output_pdf) or ".", exist_ok=True)

# ── Load pseudobulk matrices ───────────────────────────────────────────────────
Fresh    = pd.read_csv(args.fresh_csv, index_col=0)
Biopsies = pd.read_csv(args.biopsies_csv, index_col=0)

# Strip common prefix from index (e.g. "subtype_GPNMB" → "GPNMB")
Fresh.index    = Fresh.index.str.split("_").str[1]
Biopsies.index = Biopsies.index.str.split("_").str[1]

# Z-score normalise across subtypes (rows)
Fresh_scaled    = (Fresh    - Fresh.mean(axis=0))    / Fresh.std(axis=0)
Biopsies_scaled = (Biopsies - Biopsies.mean(axis=0)) / Biopsies.std(axis=0)

# Sort and align
custom_order = [
    "AIF1", "CCL3", "HIF1A", "HIST", "HSPA1A", "IFI44L", "TMEM163",
    "FRMD4A", "PICALM", "CD163", "GPNMB", "MKI67", "ERN1"
]
Fresh_scaled    = Fresh_scaled.reindex(custom_order).dropna(how="all")
Biopsies_scaled = Biopsies_scaled.reindex(custom_order).dropna(how="all")

# Shared genes
shared = Fresh_scaled.columns.intersection(Biopsies_scaled.columns)
Fresh_scaled    = Fresh_scaled[shared]
Biopsies_scaled = Biopsies_scaled[shared]

# Build full matrix: FreshMG subtypes vs ICH subtypes
corr_full  = pd.DataFrame(index=Fresh_scaled.index,
                           columns=Biopsies_scaled.index, dtype=float)
for ri in Fresh_scaled.index:
    for ci in Biopsies_scaled.index:
        r_val, _ = stats.pearsonr(Fresh_scaled.loc[ri].fillna(0),
                                   Biopsies_scaled.loc[ci].fillna(0))
        corr_full.loc[ri, ci] = r_val

# Plot
g = sns.clustermap(
    corr_full.astype(float),
    annot=False, cmap="YlOrRd",
    figsize=(6, 6),
    row_cluster=False, col_cluster=False
)
g.ax_heatmap.set_xlabel("ICH Biopsies")
g.ax_heatmap.set_ylabel("FreshMG")
g.ax_heatmap.set_title("Pearson Correlation")
g.savefig(args.output_pdf)
plt.close("all")

print(f"Saved: {args.output_pdf}")
