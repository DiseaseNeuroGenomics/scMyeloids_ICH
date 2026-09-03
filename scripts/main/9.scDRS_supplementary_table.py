#!/usr/bin/env python3
"""
Build the scDRS group-association supplementary table from the per-trait
summary CSVs written by 8.scDRS.py.

Produces the full significance table underlying Figure 2B (stroke traits) and
Supplementary Figure 3B (other brain disorders), plus BH-FDR correction, which
the figures currently do not apply.

Usage:
    python 9.scDRS_supplementary_table.py \
        --summary_dir 2024_04_01_scDRS/traits/summary \
        --out_prefix  tables/supplementary/SupplementaryTable_scDRS

Outputs:
    <out_prefix>.xlsx   sheets: README, All_traits_long, Fig2B_panel, SFig3B_panel
    <out_prefix>.csv    All_traits_long (flat, for the repo / reviewers)
"""

import argparse
import glob
import os
import re

import numpy as np
import pandas as pd

# ── Trait panels as displayed in the manuscript ───────────────────────────────
# Order here is the order the columns should appear in each figure.
FIG2B_TRAITS = [
    "ich", "ich_deep", "ich_lobar",
    "stroke", "stroke_ce", "stroke_lvd", "stroke_svd",
]

SFIG3B_TRAITS = [
    "adhd_or_asd", "alz2", "alz", "alz2noapoe", "alzKunkle", "alzWightman",
    "asd", "asd_vs_adhd", "ibd", "mdd", "mdd2", "mdd_ipsych", "ms",
    "pd", "pd_without_23andMe", "sa_mdd", "sa_mdd_bip_scz",
    "sz", "sz3", "sz_eur", "sz_mvp",
]

# Only these traits are reported. The summary directory also contains runs
# against internal / unpublished GWAS used for QC during development; those are
# not citable and are deliberately excluded from every output of this script.
REPORTED_TRAITS = FIG2B_TRAITS + SFIG3B_TRAITS

GROUP_ORDER = ["MTC_1", "MTC_2", "MTC_3", "MTC_4", "MTC_5", "MTC_6", "MTC_Prolif"]

# Human-readable trait labels for the supplementary table.
TRAIT_LABELS = {
    "ich": "Intracerebral hemorrhage (all)",
    "ich_deep": "Intracerebral hemorrhage, deep (non-lobar)",
    "ich_lobar": "Intracerebral hemorrhage, lobar",
    "stroke": "Ischemic stroke (all)",
    "stroke_ce": "Ischemic stroke, cardioembolic",
    "stroke_lvd": "Ischemic stroke, large-vessel disease",
    "stroke_svd": "Ischemic stroke, small-vessel disease",
    "alz": "Alzheimer's disease",
    "alz2": "Alzheimer's disease (alternative GWAS)",
    "alz2noapoe": "Alzheimer's disease (alternative GWAS, APOE excluded)",
    "alzKunkle": "Alzheimer's disease (Kunkle et al.)",
    "alzWightman": "Alzheimer's disease (Wightman et al.)",
    "adhd_or_asd": "ADHD or autism spectrum disorder",
    "asd": "Autism spectrum disorder",
    "asd_vs_adhd": "Autism spectrum disorder vs ADHD",
    "ibd": "Inflammatory bowel disease",
    "mdd": "Major depressive disorder",
    "mdd2": "Major depressive disorder (alternative GWAS)",
    "mdd_ipsych": "Major depressive disorder (iPSYCH)",
    "ms": "Multiple sclerosis",
    "pd": "Parkinson's disease",
    "pd_without_23andMe": "Parkinson's disease (excluding 23andMe)",
    "sa_mdd": "Suicide attempt in major depressive disorder",
    "sa_mdd_bip_scz": "Suicide attempt in MDD, bipolar disorder and schizophrenia",
    "sz": "Schizophrenia",
    "sz3": "Schizophrenia (PGC3)",
    "sz_eur": "Schizophrenia (European ancestry)",
    "sz_mvp": "Schizophrenia (Million Veteran Program)",
}

# n_ctrl is read from the summary files themselves rather than hard-coded, so
# this script stays correct if 8.scDRS.py is rerun with a different n_ctrl.
# The smallest attainable Monte-Carlo p-value is 1/(n_ctrl+1).


def load_summaries(summary_dir):
    pattern = os.path.join(summary_dir, "Metacells_identity_*_summary.csv")
    files = sorted(glob.glob(pattern))
    if not files:
        raise SystemExit(f"No summary files matched: {pattern}")

    frames = []
    skipped = []
    for f in files:
        trait = re.sub(r"^Metacells_identity_|_summary\.csv$", "",
                       os.path.basename(f))
        if trait not in REPORTED_TRAITS:
            skipped.append(trait)      # internal / unpublished GWAS - not reported
            continue
        # 8.scDRS.py writes these with sep="\t" despite the .csv extension.
        df = pd.read_csv(f, sep="\t")
        df = df.rename(columns={df.columns[0]: "group"})
        df.insert(0, "trait", trait)
        frames.append(df)

    d = pd.concat(frames, ignore_index=True)
    d["n_cell"] = d["n_cell"].astype(int)
    d["n_ctrl"] = d["n_ctrl"].astype(int)
    return d, skipped


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--summary_dir", default="2024_04_01_scDRS/traits/summary")
    ap.add_argument("--out_prefix",
                    default="tables/supplementary/SupplementaryTable_scDRS")
    args = ap.parse_args()

    d, skipped = load_summaries(args.summary_dir)

    n_ctrl_vals = sorted(d["n_ctrl"].unique())
    if len(n_ctrl_vals) > 1:
        raise SystemExit(f"Mixed n_ctrl across traits: {n_ctrl_vals}. "
                         "Rerun all traits with the same n_ctrl before tabulating.")
    n_ctrl = int(n_ctrl_vals[0])
    mcp_floor = 1.0 / (n_ctrl + 1)

    # ── Panel membership ──────────────────────────────────────────────────────
    def panel_of(t):
        return "Figure 2B" if t in FIG2B_TRAITS else "Supplementary Figure 3B"

    d["panel"] = d["trait"].map(panel_of)
    d["trait_label"] = d["trait"].map(TRAIT_LABELS).fillna(d["trait"])

    # ── Monte-Carlo p-value floor ─────────────────────────────────────────────
    # The smallest attainable assoc_mcp is 1/(n_ctrl+1). Values sitting at the
    # floor are upper bounds, not point estimates, and must be reported as
    # "< floor" rather than as exact p-values.
    floor_str = f"{mcp_floor:.3g}"
    d["assoc_mcp_at_floor"] = np.isclose(d["assoc_mcp"], mcp_floor,
                                         rtol=1e-3, atol=0.0)
    d["assoc_mcp_reported"] = np.where(
        d["assoc_mcp_at_floor"], f"<{floor_str}",
        d["assoc_mcp"].map(lambda x: f"{x:.3g}")
    )

    # ── Significance ──────────────────────────────────────────────────────────
    # Reported as the uncorrected scDRS Monte-Carlo p-value, matching the
    # asterisks in Figure 2B and Supplementary Figure 3B. See README sheet for
    # why BH correction is not applied at n_ctrl=100.
    d["significant_p05"] = d["assoc_mcp"] < 0.05

    # ── Order and write long table ────────────────────────────────────────────
    d["group"] = pd.Categorical(d["group"], categories=GROUP_ORDER, ordered=True)
    panel_rank = {"Figure 2B": 0, "Supplementary Figure 3B": 1}
    d["_p"] = d["panel"].map(panel_rank)
    order_in_panel = {t: i for i, t in enumerate(FIG2B_TRAITS + SFIG3B_TRAITS)}
    d["_t"] = d["trait"].map(order_in_panel).fillna(999)
    d = d.sort_values(["_p", "_t", "trait", "group"]).drop(columns=["_p", "_t"])

    cols = ["panel", "trait", "trait_label", "group", "n_cell", "n_ctrl",
            "assoc_mcz", "assoc_mcp", "assoc_mcp_reported", "assoc_mcp_at_floor",
            "significant_p05", "hetero_mcz", "hetero_mcp",
            "n_fdr_0.05", "n_fdr_0.1", "n_fdr_0.2"]
    cols = [c for c in cols if c in d.columns]
    d = d[cols]

    def wide(traits, value):
        sub = d[d["trait"].isin(traits)]
        w = sub.pivot_table(index="group", columns="trait", values=value,
                            observed=False)
        return w.reindex(index=GROUP_ORDER, columns=traits)

    os.makedirs(os.path.dirname(args.out_prefix) or ".", exist_ok=True)
    d.to_csv(f"{args.out_prefix}.csv", index=False)

    # The justification for not correcting depends on the Monte-Carlo resolution:
    # a BH q-value below 0.05 across k tests needs at least k/5 tests at the floor.
    largest_panel = int(d.groupby("panel").size().max())
    if mcp_floor * largest_panel / max(1, int(np.ceil(largest_panel / 5))) > 0.05:
        mtc_text = (
            "None. Significance is the uncorrected scDRS Monte-Carlo p-value "
            "(assoc_mcp < 0.05), as marked in Figure 2B and Supplementary Figure 3B. "
            f"Benjamini-Hochberg correction is uninformative at n_ctrl={n_ctrl}: the "
            f"p-value floor is {floor_str}, so across the largest panel "
            f"({largest_panel} tests) a q-value below 0.05 would require at least "
            f"{int(np.ceil(largest_panel / 5))} tests to sit simultaneously at the "
            "floor. The procedure would therefore be bounded by the resolution of the "
            "permutation rather than by the data. scDRS group-level associations are "
            "reported as exploratory and interpreted together with the effect "
            "direction (assoc_mcz), not as confirmatory tests.")
    else:
        mtc_text = (
            "None applied in this table. Significance is the uncorrected scDRS "
            "Monte-Carlo p-value (assoc_mcp < 0.05), matching the figure markers. "
            f"At n_ctrl={n_ctrl} the p-value floor is {floor_str}, so Benjamini-Hochberg "
            "correction across a panel is now feasible and should be considered if the "
            "associations are to be presented as confirmatory rather than exploratory.")

    readme = pd.DataFrame({"Field": [
        "Source", "Method", "n_ctrl", "assoc_mcz", "assoc_mcp",
        "Monte-Carlo p floor", "Multiple-testing correction",
        "hetero_mcp", "Figure 2B traits", "Supplementary Figure 3B traits",
        "Scope",
    ], "Description": [
        "scDRS group-level association test (scdrs.method.downstream_group_analysis), "
        "grouping variable Metacells_identity; generated by 8.scDRS.py, seed 12345678.",
        f"scdrs.score_cell(ctrl_match_key='mean_var', n_ctrl={n_ctrl}, "
        "weight_opt='vs'); gene sets from MAGMA top-1000 genes per trait.",
        str(n_ctrl),
        "Monte-Carlo z-score for the group-level association; positive values "
        "indicate higher disease-relevance scores in that metacell cluster.",
        "Monte-Carlo p-value for the group-level association (uncorrected).",
        f"With n_ctrl={n_ctrl} the smallest attainable assoc_mcp is 1/(n_ctrl+1) = "
        f"{floor_str}. Rows flagged assoc_mcp_at_floor=TRUE are upper bounds and "
        f"are reported as '<{floor_str}', not as exact p-values.",
        mtc_text,
        "Monte-Carlo p-value for within-group heterogeneity of disease-relevance "
        "scores. Reported for completeness; not used to define significance.",
        ", ".join(FIG2B_TRAITS),
        ", ".join(SFIG3B_TRAITS),
        "This table reports only the traits displayed in Figure 2B and Supplementary "
        "Figure 3B. Additional gene sets were scored during method development for "
        "internal quality control; those are not included here.",
    ]})

    with pd.ExcelWriter(f"{args.out_prefix}.xlsx", engine="openpyxl") as xl:
        readme.to_excel(xl, sheet_name="README", index=False)
        d.to_excel(xl, sheet_name="All_traits_long", index=False)
        wide(FIG2B_TRAITS, "assoc_mcz").to_excel(xl, sheet_name="Fig2B_zscore")
        wide(FIG2B_TRAITS, "assoc_mcp").to_excel(xl, sheet_name="Fig2B_p")
        wide(SFIG3B_TRAITS, "assoc_mcz").to_excel(xl, sheet_name="SFig3B_zscore")
        wide(SFIG3B_TRAITS, "assoc_mcp").to_excel(xl, sheet_name="SFig3B_p")

    # ── Console report ────────────────────────────────────────────────────────
    print(f"Reported traits: {d['trait'].nunique()}  |  tests: {len(d)}")
    if skipped:
        print(f"Excluded {len(skipped)} non-reported gene set(s) present in "
              f"{args.summary_dir} (internal QC runs).")

    for panel in ["Figure 2B", "Supplementary Figure 3B"]:
        sub = d[d["panel"] == panel]
        if not len(sub):
            continue
        print(f"\n{panel}: {sub['trait'].nunique()} traits, {len(sub)} tests")
        print(f"  significant (assoc_mcp < 0.05) : "
              f"{int(sub['significant_p05'].sum())}")
        print(f"  at Monte-Carlo p floor         : "
              f"{int(sub['assoc_mcp_at_floor'].sum())} "
              f"(report as P < {floor_str})")
        for _, r in sub[sub["significant_p05"]].iterrows():
            print(f"    {r['group']:<11} {r['trait']:<20} "
                  f"z={r['assoc_mcz']:+.2f}  P={r['assoc_mcp_reported']}")

    missing = [t for t in REPORTED_TRAITS if t not in set(d["trait"])]
    if missing:
        print(f"\nWARNING - expected but not found on disk: {missing}")

    print(f"\nWrote {args.out_prefix}.xlsx and {args.out_prefix}.csv")


if __name__ == "__main__":
    main()
