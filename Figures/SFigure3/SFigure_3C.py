#!/usr/bin/env python3
# Author: Dimitrios Kyriakis
# Description: scDRS enrichment heatmap for non-stroke neuropsychiatric traits
#              across myeloid metacell clusters. Produces SFigure 3C.
# Usage: python SFigure_3C.py
#   --summary_dir  Directory containing Metacells_identity_<trait>_summary.csv files
#   --output_pdf   Output PDF path (default: SFigures/SFigure3/2024_12_Sup.Fig_3C.pdf)

import argparse
import os

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser(
    description="scDRS non-stroke trait enrichment heatmap"
)
parser.add_argument("--summary_dir", required=True,
                    help="Directory with *_summary.csv files")
parser.add_argument("--output_pdf",  default="SFigures/SFigure3/2024_12_Sup.Fig_3C.pdf",
                    help="Output PDF path")
args = parser.parse_args()

os.makedirs(os.path.dirname(args.output_pdf) or ".", exist_ok=True)

# ── Trait lists (non-stroke) ──────────────────────────────────────────────────
alz_columns = ["alz2", "alz2noapoe", "alzWightman", "alz", "alzKunkle"]
scz_columns = ["sz3", "sz_mvp", "sz", "sa_mdd_bip_scz", "sz_eur"]
mdd_columns = ["sa_mdd", "mdd", "mdd2", "mdd_ipsych"]
asd_columns = ["adhd_or_asd", "asd", "asd_vs_adhd"]
pd_columns  = ["pd", "pd_without_23andMe"]
ms_columns  = ["ms"]
bd_columns  = ["ibd"]

not_stroke = (alz_columns + scz_columns + mdd_columns +
              asd_columns + pd_columns + ms_columns + bd_columns)

tbls_z, tbls_p, names_ok = [], [], []
for trait in not_stroke:
    path = os.path.join(args.summary_dir,
                        f"Metacells_identity_{trait}_summary.csv")
    if not os.path.exists(path):
        print(f"WARNING: missing {path} — skipping")
        continue
    tbl = pd.read_table(path, sep="\t", on_bad_lines="skip")
    tbls_z.append(tbl.set_index("group")["assoc_mcz"].rename(trait))
    tbls_p.append(tbl.set_index("group")["assoc_mcp"].rename(trait))
    names_ok.append(trait)

if not tbls_z:
    raise SystemExit("No trait summary files found — check --summary_dir")

df_z = pd.concat(tbls_z, axis=1).sort_index(axis=1)
df_p = pd.concat(tbls_p, axis=1).sort_index(axis=1)

annot = df_p.apply(lambda x: ["*" if y <= 0.05 else "" for y in x])

sns.set(font_scale=1.2)
ax = sns.clustermap(
    df_z, metric="euclidean",
    square=True, figsize=(15, 5),
    col_cluster=False, row_cluster=False,
    yticklabels=True, cmap="RdBu_r",
    annot=annot, fmt="",
    cbar_pos=(1, 0.55, 0.03, 0.2)
)
plt.setp(ax.ax_heatmap.get_xticklabels(), rotation=90)
ax.savefig(args.output_pdf)
plt.savefig(args.output_pdf.replace(".pdf", ".png"))
plt.close("all")

print(f"Saved: {args.output_pdf}")
