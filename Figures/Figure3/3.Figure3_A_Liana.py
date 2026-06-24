#!/usr/bin/env python3
# Author: Dimitrios Kyriakis
# Description: LIANA cell-cell interaction dot plots for Figure 3A-B.
#              Loads pre-computed LIANA results from h5ad and differential CCI CSV,
#              filters to |log2FC| > 1.5, and plots source/target panels.
# Usage: python 3.Figure3_A_Liana.py
#   --h5ad        Path to 2025_Figure3B_Liana.h5ad
#   --cci_csv     Path to 2024_03_28_Log2FC_CCIs_Targets_MTC123_vs_MTC456.csv
#   --output_dir  Output directory (default: current dir)
#   --lfc_thresh  |log2FC| threshold for filtering (default: 1.5)

import argparse
import os

import numpy as np
import pandas as pd
import scanpy as sc
import plotnine as p9
from plotnine import (ggplot, geom_point, aes, facet_grid, facet_wrap,
                      labs, theme_bw, theme, element_text, element_rect,
                      scale_size_continuous, scale_color_gradient2,
                      scale_color_cmap, coord_flip, theme_minimal)

parser = argparse.ArgumentParser(description="LIANA CCI dot plots")
parser.add_argument("--h5ad",       required=True,
                    help="Path to 2025_Figure3B_Liana.h5ad")
parser.add_argument("--cci_csv",    required=True,
                    help="Path to differential CCI CSV")
parser.add_argument("--output_dir", default=".",
                    help="Output directory")
parser.add_argument("--lfc_thresh", type=float, default=1.5,
                    help="|log2FC| threshold (default: 1.5)")
args = parser.parse_args()

os.makedirs(args.output_dir, exist_ok=True)

# ── Load data ─────────────────────────────────────────────────────────────────
adata   = sc.read_h5ad(args.h5ad)
cci_df  = pd.read_csv(args.cci_csv)

# Filter to differential interactions
sig_cci = cci_df[cci_df["log2FC"].abs() > args.lfc_thresh].copy()

liana_res = adata.uns["liana_res"].copy()
liana_res["interaction"] = (liana_res["ligand_complex"] + "."
                             + liana_res["receptor_complex"])

sig_interactions = set(sig_cci["CCI"].str.split(":").str[1]
                        if "CCI" in sig_cci.columns else
                        sig_cci["ligand_complex"] + "." + sig_cci["receptor_complex"])

liana_filt = liana_res[liana_res["interaction"].isin(sig_interactions)].copy()
liana_filt["neg_mag"] = 1 - liana_filt["magnitude_rank"]
liana_filt["neg_spc"] = 1 - liana_filt["specificity_rank"]


def dotplot4(df, source_label, orientation="source", figure_size=(10, 8),
             output_path=None):
    """Dot plot coloured by magnitude, sized by specificity rank."""
    sub = df[df["source" if orientation == "source" else "target"] == source_label]
    if sub.empty:
        return
    p = (
        ggplot(sub, aes(x="target" if orientation == "source" else "source",
                        y="interaction",
                        colour="neg_mag", size="neg_spc"))
        + geom_point()
        + scale_color_cmap("PuRd")
        + scale_size_continuous(range=(1, 8))
        + facet_wrap("~source" if orientation == "source" else "~target")
        + labs(x="Target cell", y="Interaction (Ligand.Receptor)",
               colour="Magnitude\n(1-rank)", size="Specificity\n(1-rank)")
        + theme_bw()
        + theme(
            axis_text_x=element_text(size=10, angle=90),
            axis_text_y=element_text(size=9),
            strip_background=element_rect(fill="white"),
            strip_text=element_text(size=12),
            figure_size=figure_size
        )
    )
    if output_path:
        p.save(output_path, dpi=150)
        print(f"Saved: {output_path}")
    return p


# Source panels (MTC_123 and MTC_456 as sources)
for mtc in ["MTC_123", "MTC_456"]:
    out = os.path.join(args.output_dir, f"Figure3B_{mtc}_source.pdf")
    dotplot4(liana_filt, mtc, orientation="source", output_path=out)

# Target panels
for mtc in ["MTC_123", "MTC_456"]:
    out = os.path.join(args.output_dir, f"Figure3B_{mtc}_target.pdf")
    dotplot4(liana_filt, mtc, orientation="target", output_path=out)
