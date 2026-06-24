#!/usr/bin/env python3
# Author: Dimitrios Kyriakis
# Description: LIANA cell-cell interaction dot plot for SFigure 5.
#              Filters to IL/IFN/complement interactions and plots dot plot
#              coloured by magnitude rank, sized by specificity rank.
# Usage: python SFigure5.py
#   --h5ad        Path to immune metacell h5ad with liana_res in uns
#   --output_pdf  Output PDF path (default: SFigures/SFigure5/SFigure5.pdf)

import argparse
import os

import numpy as np
import pandas as pd
import scanpy as sc
import plotnine as p9
from plotnine import (ggplot, geom_point, aes, facet_grid, labs,
                      theme_bw, theme, element_text, element_rect,
                      scale_size_continuous, scale_color_cmap)

parser = argparse.ArgumentParser(description="LIANA dotplot for SFigure 5")
parser.add_argument("--h5ad",       required=True,
                    help="Path to h5ad with liana_res in uns")
parser.add_argument("--output_pdf", default="SFigures/SFigure5/SFigure5.pdf",
                    help="Output PDF path")
args = parser.parse_args()

os.makedirs(os.path.dirname(args.output_pdf) or ".", exist_ok=True)

adata   = sc.read_h5ad(args.h5ad)
CCI_Df  = adata.uns["liana_res"].copy()

# Filter to IL/IFN/complement interactions
start_with_cci = ("IL10", "IFN", "IL", "C3", "C5")
mask = (CCI_Df["receptor_complex"].str.startswith(start_with_cci) |
        CCI_Df["ligand_complex"].str.startswith(start_with_cci))
subset_df = CCI_Df[mask].copy()
subset_df["interaction"] = (subset_df["ligand_complex"] + "."
                             + subset_df["receptor_complex"])
subset_df["neg_mag"] = 1 - subset_df["magnitude_rank"]
subset_df["neg_spc"] = 1 - subset_df["specificity_rank"]

# Retain only significant interactions
subset_df = subset_df[subset_df["cellphone_pvals"] < 0.05]

if subset_df.empty:
    print("No significant interactions after filtering — nothing to plot.")
else:
    p = (
        ggplot(subset_df,
               aes(x="target", y="interaction",
                   colour="neg_mag", size="neg_spc"))
        + geom_point()
        + facet_grid("~source")
        + scale_size_continuous(range=(1, 8))
        + scale_color_cmap("PuRd")
        + labs(
            x      = "Target cell",
            y      = "Interactions (Ligand.Receptor)",
            colour = "Magnitude\n(1-rank)",
            size   = "Specificity\n(1-rank)",
            title  = "Source"
        )
        + theme_bw()
        + theme(
            legend_text      = element_text(size=12),
            strip_background = element_rect(fill="white"),
            strip_text       = element_text(size=13, colour="black"),
            axis_text_y      = element_text(size=9, colour="black"),
            axis_text_x      = element_text(size=10, angle=90),
            figure_size      = (14, 8)
        )
    )
    p.save(args.output_pdf, dpi=150)
    print(f"Saved: {args.output_pdf}")
