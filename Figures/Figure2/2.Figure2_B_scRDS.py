#!/usr/bin/env python3
# Author: Dimitrios Kyriakis
# Description: scDRS enrichment heatmaps for stroke traits (Figure 2B)
#              and non-stroke neuropsychiatric traits (SFigure 3C).
#              Reads pre-computed per-trait summary CSVs.
# Usage: python 2.Figure2_B_scRDS.py
#   --summary_dir   Directory containing *_summary.csv files
#   --output_fig    Output PDF for stroke heatmap (default: Figure2B.pdf)
#   --output_sfig   Output PDF for non-stroke heatmap (default: SFigure3C.pdf)

import argparse
import os
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser(
    description="scDRS enrichment heatmaps from summary CSVs"
)
parser.add_argument("--summary_dir", required=True,
                    help="Directory with Metacells_identity_<trait>_summary.csv files")
parser.add_argument("--output_fig",  default="Figure2B.pdf",
                    help="Output PDF for stroke traits heatmap")
parser.add_argument("--output_sfig", default="SFigure3C.pdf",
                    help="Output PDF for non-stroke traits heatmap")
args = parser.parse_args()

# ── Trait lists ───────────────────────────────────────────────────────────────
stroke_columns = ["ich", "ich_deep", "ich_lobar",
                  "stroke", "stroke_ce", "stroke_lvd", "stroke_svd"]

alz_columns = ["alz2", "alz2noapoe", "alzWightman", "alz", "alzKunkle"]
scz_columns = ["sz3", "sz_mvp", "sz", "sa_mdd_bip_scz", "sz_eur"]
mdd_columns = ["sa_mdd", "mdd", "mdd2", "mdd_ipsych"]
asd_columns = ["adhd_or_asd", "asd", "asd_vs_adhd"]
pd_columns  = ["pd", "pd_without_23andMe"]
ms_columns  = ["ms"]
bd_columns  = ["ibd"]

not_stroke = (alz_columns + scz_columns + mdd_columns +
              asd_columns + pd_columns + ms_columns + bd_columns)


def _load_traits(trait_list, summary_dir):
    """Return (zscore_df, pval_df) for a list of traits."""
    tbls_z, tbls_p, nms = [], [], []
    for trait in trait_list:
        path = os.path.join(summary_dir,
                            f"Metacells_identity_{trait}_summary.csv")
        if not os.path.exists(path):
            print(f"WARNING: missing {path} — skipping")
            continue
        tbl = pd.read_table(path, sep="\t", on_bad_lines="skip")
        tbls_z.append(tbl["assoc_mcz"])
        tbls_p.append(tbl["assoc_mcp"])
        nms.append(trait)
    if not tbls_z:
        return None, None
    first_tbl = pd.read_table(
        os.path.join(summary_dir, f"Metacells_identity_{nms[0]}_summary.csv"),
        sep="\t", on_bad_lines="skip")
    groups = first_tbl["group"].values
    df_z = pd.DataFrame({n: pd.read_table(
        os.path.join(summary_dir, f"Metacells_identity_{n}_summary.csv"),
        sep="\t", on_bad_lines="skip")["assoc_mcz"].values
        for n in nms}, index=groups)
    df_p = pd.DataFrame({n: pd.read_table(
        os.path.join(summary_dir, f"Metacells_identity_{n}_summary.csv"),
        sep="\t", on_bad_lines="skip")["assoc_mcp"].values
        for n in nms}, index=groups)
    return df_z, df_p


def _plot_heatmap(df_z, df_p, output_path, figsize, col_cluster=False, row_cluster=False):
    annot = df_p.apply(lambda x: ["*" if y <= 0.05 else "" for y in x])
    sns.set(font_scale=1.2)
    ax = sns.clustermap(
        df_z, metric="euclidean", square=True, figsize=figsize,
        col_cluster=col_cluster, row_cluster=row_cluster,
        yticklabels=True, cmap="RdBu_r",
        annot=annot, fmt="",
        cbar_pos=(1, 0.55, 0.03, 0.2)
    )
    plt.setp(ax.ax_heatmap.get_xticklabels(), rotation=90)
    ax.savefig(output_path)
    plt.close("all")
    print(f"Saved: {output_path}")


# ── Figure 2B — stroke traits ─────────────────────────────────────────────────
df_z_stroke, df_p_stroke = _load_traits(stroke_columns, args.summary_dir)
if df_z_stroke is not None:
    _plot_heatmap(df_z_stroke, df_p_stroke, args.output_fig,
                  figsize=(4, 4), col_cluster=False, row_cluster=False)

# ── SFigure 3C — non-stroke traits ───────────────────────────────────────────
df_z_ns, df_p_ns = _load_traits(not_stroke, args.summary_dir)
if df_z_ns is not None:
    df_z_ns = df_z_ns.sort_index(axis=1)
    df_p_ns = df_p_ns.sort_index(axis=1)
    _plot_heatmap(df_z_ns, df_p_ns, args.output_sfig,
                  figsize=(15, 5), col_cluster=False, row_cluster=False)
