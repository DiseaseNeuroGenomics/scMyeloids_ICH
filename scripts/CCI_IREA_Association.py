#!/usr/bin/env python3
# Author: Dimitrios Kyriakis
# Description: Fisher exact test for overlap between significant cell-cell
#              interactions (LIANA) and immune-response enriched cytokines (IREA)
#              in MTC_123 vs MTC_456 metaclusters.
# Usage: python CCI_IREA_Association.py
#   --h5ad           Path to h5ad with liana_res in uns
#   --irea_123       Path to MTC_123 IREA CSV (columns: Cytokine, padj)
#   --irea_456       Path to MTC_456 IREA CSV (columns: Cytokine, padj)
#   --output_csv     Output results CSV (default: CCI_IREA_overlap.csv)
#   --padj_thresh    Adjusted p-value threshold for both sets (default: 0.05)

import argparse
import os

import numpy as np
import pandas as pd
import scipy.stats as stats

parser = argparse.ArgumentParser(
    description="Fisher exact test: LIANA CCI vs IREA cytokine overlap"
)
parser.add_argument("--h5ad",       required=True,
                    help="h5ad file with liana_res in uns")
parser.add_argument("--irea_123",   required=True,
                    help="MTC_123 IREA CSV")
parser.add_argument("--irea_456",   required=True,
                    help="MTC_456 IREA CSV")
parser.add_argument("--output_csv", default="CCI_IREA_overlap.csv",
                    help="Output results CSV")
parser.add_argument("--padj_thresh", type=float, default=0.05,
                    help="Adjusted p-value threshold (default: 0.05)")
args = parser.parse_args()

import scanpy as sc
adata   = sc.read_h5ad(args.h5ad)
cci_df  = adata.uns["liana_res"].copy()

irea_123 = pd.read_csv(args.irea_123)
irea_456 = pd.read_csv(args.irea_456)
irea_123["Cytokine"] = irea_123["Cytokine"].str.upper()
irea_456["Cytokine"] = irea_456["Cytokine"].str.upper()

# ── Significant sets ──────────────────────────────────────────────────────────
mtc_123_cci = cci_df[
    (cci_df["source"] == "MTC_123") &
    (cci_df["cellphone_pvals"] < 0.05)
]
mtc_456_cci = cci_df[
    (cci_df["source"] == "MTC_456") &
    (cci_df["cellphone_pvals"] < 0.05)
]

mtc_123_cci_genes = (set(mtc_123_cci["ligand_complex"])
                     | set(mtc_123_cci["receptor_complex"]))
mtc_456_cci_genes = (set(mtc_456_cci["ligand_complex"])
                     | set(mtc_456_cci["receptor_complex"]))

mtc_123_irea_sig = set(
    irea_123.query(f"padj < {args.padj_thresh}")["Cytokine"]
)
mtc_456_irea_sig = set(
    irea_456.query(f"padj < {args.padj_thresh}")["Cytokine"]
)

background = (set(cci_df["ligand_complex"])
              | set(cci_df["receptor_complex"])
              | set(irea_456["Cytokine"])
              | set(irea_123["Cytokine"]))
N = len(background)

results = []
for label, cci_genes, irea_genes in [
    ("MTC_123", mtc_123_cci_genes, mtc_123_irea_sig),
    ("MTC_456", mtc_456_cci_genes, mtc_456_irea_sig),
]:
    a = len(cci_genes & irea_genes)
    b = len(cci_genes - irea_genes)
    c = len(irea_genes - cci_genes)
    d = N - (a + b + c)
    _, pval = stats.fisher_exact([[a, b], [c, d]])
    results.append({
        "Cluster":           label,
        "N_CCI_sig":         len(cci_genes),
        "N_IREA_sig":        len(irea_genes),
        "Overlap":           a,
        "Overlap_genes":     "; ".join(sorted(cci_genes & irea_genes)),
        "Fisher_p":          pval,
        "Background_size":   N,
    })
    print(f"{label}: CCI={len(cci_genes)}, IREA={len(irea_genes)}, "
          f"Overlap={a}, p={pval:.4g}")

res_df = pd.DataFrame(results)
res_df.to_csv(args.output_csv, index=False)
print(f"Saved: {args.output_csv}")
