# =============================================================================
# Figure 3B — LIANA cell-cell interactions
# Writes: Table16_CCI_LIANA.csv  and  Figure3_B_2026_08_12.pdf
# Run from inside the Liana/ folder.
# =============================================================================
import sys
import logging
import numpy as np
import scanpy as sc
import hdf5plugin
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle

plt.rcParams["pdf.fonttype"] = 42
plt.rcParams["ps.fonttype"] = 42

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s | %(levelname)s | %(message)s",
    handlers=[logging.StreamHandler(sys.stderr)]
)
logger = logging.getLogger(__name__)
H5AD      = "2024_04_22_Immune_Metacells_combined_stored_counts.h5ad.ztsd"
TABLE_OUT = "Table16_CCI_LIANA.csv"
FIG_OUT   = "Figure3_B_2026_08_12.pdf"
SPEC_CUT, MAG_CUT = 0.05, 0.05
MIN_DOT_AREA = 5      
MAX_DOT_AREA = 150    
SOURCES = ["MTC_123", "MTC_456", "MTC_Prolif", "MONO", "NEUT"]
TARGETS = ["MTC_123", "MTC_456", "MTC_Prolif", "MONO", "NEUT", "T Cells", "B Cells"]
PAIRS = [("APOE", "LRP1"), ("APOE", "SDC2"), ("APOE", "SORL1"), ("APOE", "TREM2"),
         ("C3", "C3AR1"), ("C3", "CD81"), ("C3", "IFITM1"), ("C3", "ITGAX"), ("C3", "LRP1"),
         ("CD14", "ITGB1"), ("CD14", "ITGB2"), ("CD14", "TLR4"),
         ("IL10", "IL10RA_IL10RB"),
         ("LGALS1", "CD69"), ("LGALS1", "ITGB1"),
         ("SPP1", "CD44"), ("SPP1", "ITGA4_ITGB1"), ("SPP1", "ITGAV_ITGB1"),
         ("SPP1", "ITGAV_ITGB5"), ("SPP1", "PTGER4")]
COL = {"MTC_123": "#7F7F7F", "MTC_456": "#8B1A1A", "MTC_Prolif": "#F0857D",
       "MONO": "#7FB2DC", "NEUT": "#DCDC8C", "T Cells": "#8CC63F", "B Cells": "#F5A83C"}
def scale_area(val, val_min, val_max, area_min=MIN_DOT_AREA, area_max=MAX_DOT_AREA):
    if val_max == val_min:
        return np.full_like(val, area_max)
    val_clipped = np.clip(val, val_min, val_max)
    norm = (val_clipped - val_min) / (val_max - val_min)
    return norm * (area_max - area_min) + area_min
# ----------------------------------------------------------- 1. the table ----
logger.info("Reading AnnData...")
adata = sc.read_h5ad(H5AD)
table = adata.uns["liana_res"].copy()
table = table[(table.specificity_rank <= SPEC_CUT) &
              (table.magnitude_rank  <= MAG_CUT)].sort_values("magnitude_rank")
table.to_csv(TABLE_OUT, index=False)
logger.info(f"Saved {TABLE_OUT}: {len(table)} interactions")
# ------------------------------------------------------- 2. the plot data ----
df = table.copy()
df["interaction"] = df.ligand_complex + " -> " + df.receptor_complex
y_order = [f"{l} -> {r}" for l, r in PAIRS]
ypos    = {k: i for i, k in enumerate(y_order)}
df = df[df.interaction.isin(y_order) &
        df.source.isin(SOURCES) &
        df.target.isin(TARGETS)].copy()
df["spec"] = -np.log10(df.specificity_rank.clip(lower=1e-13))
df["mag"]  = -np.log10(df.magnitude_rank.clip(lower=1e-13))
min_spec = df["spec"].min()
max_spec = df["spec"].max()
max_mag  = df["mag"].max()
logger.info(f"Plotting {len(df)} points.")
cmap = LinearSegmentedColormap.from_list(
    "mag", ["#FFFFFF", "#F2C6DE", "#D954A0", "#8C1A5B", "#3B0A28"])
norm = Normalize(0, max_mag)
# ----------------------------------------------------------- 3. the panel ----
NY = len(y_order)
# Widened to 14, shortened to 4.8
fig, axes = plt.subplots(1, len(SOURCES), figsize=(14, 4.8), sharey=True,
                         gridspec_kw={"wspace": 0.10})
# Adjusted margins to leave room for the outside bars
fig.subplots_adjust(left=0.20, right=0.72, top=0.86, bottom=0.18)
for ax, src in zip(axes, SOURCES):
    sub = df[df.source == src]
    
    dot_areas = scale_area(sub.spec, val_min=min_spec, val_max=max_spec)
    
    ax.scatter([TARGETS.index(t) for t in sub.target],
               [ypos[i] for i in sub.interaction],
               s=dot_areas, c=sub.mag, cmap=cmap, norm=norm,
               edgecolors="none", zorder=3)
    # Tightened limits to only enclose the grid/dots
    ax.set_xlim(-0.5, len(TARGETS) - 0.5)
    ax.set_ylim(-0.5, NY - 0.5)
    ax.set_xticks([])
    ax.grid(True, color="#EFEFEF", lw=0.5)
    ax.set_axisbelow(True)
    
    for s_spine in ax.spines.values():
        s_spine.set_color("#BBBBBB")
    # Drawn OUTSIDE the axes (clip_on=False)
    # Source bar (Top)
    ax.add_patch(Rectangle((-0.5, NY - 0.5 + 0.2), len(TARGETS), 0.6,
                           color=COL[src], clip_on=False, zorder=4))
    # Target bars (Bottom)
    for j, t in enumerate(TARGETS):
        ax.add_patch(Rectangle((j - 0.4, -0.5 - 0.8), 0.8, 0.6,
                               color=COL[t], clip_on=False, zorder=4))
axes[0].set_yticks(range(NY))
axes[0].set_yticklabels(y_order, fontsize=8)
axes[0].set_ylabel("Interactions (Ligand -> Receptor)", fontsize=9)
# Shifted text slightly to accommodate the new layout
fig.text(0.46, 0.95, "Source",  ha="center", fontsize=10)
fig.text(0.46, 0.05, "Targets", ha="center", fontsize=10)
# ---------------------------------------------------------- 4. the legends ---
legend_breaks = np.linspace(min_spec, max_spec, 5)
size_handles = [
    Line2D([], [], marker="o", ls="", color="black",
           markersize=np.sqrt(scale_area(val, val_min=min_spec, val_max=max_spec)), 
           label=f"{val:.1f}")
    for val in legend_breaks
]
leg = fig.legend(handles=size_handles, title="Specificity rank\n(-log10)",
                 loc="upper left", bbox_to_anchor=(0.765, 0.90),
                 frameon=False, fontsize=8, title_fontsize=8, labelspacing=1.3)
fig.add_artist(leg)
fig.legend(handles=[Rectangle((0, 0), 1, 1, color=COL[k], label=k) for k in TARGETS],
           loc="upper left", bbox_to_anchor=(0.765, 0.45),
           frameon=False, fontsize=8)
cax = fig.add_axes([0.91, 0.18, 0.015, 0.20])
cb = fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap), cax=cax)
cb.set_label("Magnitude rank (-log10)", fontsize=8)
cb.ax.tick_params(labelsize=7)
plt.savefig(FIG_OUT)
logger.info(f"Successfully wrote figure to {FIG_OUT}")
