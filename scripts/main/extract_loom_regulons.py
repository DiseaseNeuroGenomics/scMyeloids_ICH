#!/usr/bin/env python3
# ============================================================================
# extract_loom_regulons.py
# Project: Decoding Immune Landscapes and Therapeutic Targets in Intracerebral Hemorrhage
# Author:  Dimitrios Kyriakis
# ============================================================================

"""
Extract the RegulonsAUC (cells x regulons) and Regulons (genes x regulons,
binary incidence) attributes from a pySCENIC AUCell-output loom, and dump
them to plain CSV.

Workaround for R's SCopeLoomR failing to read these specific compound-typed
HDF5 attributes from a loom built via the current SCopeLoomR GitHub `master`
(get_regulons_AUC() reports "attribute not available" even though `h5ls`
confirms the data is physically present) — the file that wrote it (Python
loompy/pyscenic) can read its own format reliably, so this bypasses R's
reader entirely for this one file rather than chasing the exact HDF5
compound-type interop mismatch further.

Usage (must run inside the pyscenic singularity container, since loompy/
pandas there are what actually wrote this file):
  singularity run --bind <FIGURE_CODES_DIR>:<FIGURE_CODES_DIR> \
    <FIGURE_CODES_DIR>/database/aertslab-pyscenic-0.12.1.sif \
    python <FIGURE_CODES_DIR>/pipelines/extract_loom_regulons.py \
      <loom_path> <output_prefix>

Output:
  <output_prefix>_RegulonsAUC.csv  (rows = cells, columns = regulons)
  <output_prefix>_Regulons.csv     (rows = genes, columns = regulons, 0/1)
"""

import sys

import loompy as lp
import pandas as pd


def main():
    if len(sys.argv) != 3:
        print("Usage: extract_loom_regulons.py <loom_path> <output_prefix>", file=sys.stderr)
        sys.exit(1)

    loom_path, output_prefix = sys.argv[1], sys.argv[2]

    with lp.connect(loom_path, mode="r", validate=False) as ds:
        auc_mtx = pd.DataFrame(ds.ca["RegulonsAUC"], index=ds.ca["CellID"])
        auc_out = f"{output_prefix}_RegulonsAUC.csv"
        auc_mtx.to_csv(auc_out)
        print(f"Saved {auc_out}: {auc_mtx.shape[0]} cells x {auc_mtx.shape[1]} regulons")

        regulons_mtx = pd.DataFrame(ds.ra["Regulons"], index=ds.ra["Gene"])
        reg_out = f"{output_prefix}_Regulons.csv"
        regulons_mtx.to_csv(reg_out)
        print(f"Saved {reg_out}: {regulons_mtx.shape[0]} genes x {regulons_mtx.shape[1]} regulons")


if __name__ == "__main__":
    main()
