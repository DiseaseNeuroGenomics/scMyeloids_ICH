# Decoding Immune Landscapes and Therapeutic Targets in Intracerebral Hemorrhage

![Project Overview](Overview.jpeg)

> Kyriakis D, Wang X, Pavlopoulos A, Vicari J, Kleopoulos SP, Shao Z, Fullard JF,
> Lee D, Georgakopoulos A, Carvalho Poyraz F, Voloudakis G, Skupin A, Kellner CP,
> Roussos P. *Journal TBD.*
>
> Corresponding author: Panos Roussos (panagiotis.roussos@mssm.edu)

---

## Abstract

Intracerebral hemorrhage (ICH) is the most fatal type of stroke. After ICH,
mitigating secondary brain injury driven by the immune response is the key to
improving patient outcomes and reducing mortality rates, and yet, we don't have
effective treatments. Using single-cell RNA sequencing of myeloid cells from
30 patients with acute ICH, we identified transcription factors and distinct
myeloid activation states associated with favorable or unfavorable outcomes at
6 months post-ICH. Unfavorable outcomes were associated with genes such as
PPARγ and SPP1, as well as pathways related to lipid metabolism and
neuroinflammation, whereas protective states were characterized by complement
pathway activity. Cell-cell interaction analysis highlighted key ligand-receptor
signaling differences, including SPP1→ITGA4/ITGB1 and APOE→TREM2/SORL1
interactions linked to unfavorable outcomes. Transcriptional networks underlying
immune cell polarization revealed dynamic shifts in myeloid activation states
that correlated with recovery or deterioration after ICH. Computational drug
repurposing identified all-trans retinoic acid (RA) and Everolimus as promising
therapeutic candidates. Consistently, in an in vivo ICH model, treatment with
Everolimus or RA led to significant improvements in motor function recovery, as
demonstrated by enhanced rotarod performance at all measured time points. These
findings offer critical insights into the immune pathways shaping ICH outcomes
and establish Everolimus and retinoic acid as promising therapeutic interventions
in a preclinical ICH model.

---

## Authors

| Author | Affiliations |
|--------|-------------|
| Dimitrios Kyriakis | 1, 2, 3, 4, 5 |
| Xinyi Wang | 1, 3, 4, 5 |
| Angelos Pavlopoulos | — |
| James Vicari | 1, 3, 4, 5 |
| Steve P. Kleopoulos | 1, 3, 4, 5 |
| Zhiping Shao | 1, 3, 4, 5 |
| John F. Fullard | 1, 3, 4, 5 |
| Donghoon Lee | 1, 3, 4, 5 |
| Anastasios Georgakopoulos | — |
| Fernanda Carvalho Poyraz | — |
| Georgios Voloudakis | — |
| Alex Skupin | 2 |
| Christopher P. Kellner | — |
| Panos Roussos | 1, 3, 4, 5, 6 |

**Affiliations**
1. Friedman Brain Institute, Icahn School of Medicine at Mount Sinai, New York, NY, USA
2. Luxembourg Center for Systems Biomedicine (LCSB), University of Luxembourg, Esch-sur-Alzette, Luxembourg
3. Center for Disease Neurogenomics, Icahn School of Medicine at Mount Sinai, New York, NY, USA
4. Department of Psychiatry, Icahn School of Medicine at Mount Sinai, New York, NY, USA
5. Department of Genetics and Genomic Sciences, Icahn School of Medicine at Mount Sinai, New York, NY, USA
6. Mental Illness Research, Education and Clinical Centers, James J. Peters VA Medical Center, Bronx, New York, USA

---

## Interactive data browser

Single-cell data with UMAP embeddings, cell-type annotations, and per-cell
metadata are available via CellXGene:
<https://cellxgene.cziscience.com/e/a34a4892-8ec2-4330-9198-81fc31d034e5.cxg/>

---

## Repository structure

```
├── scripts/                     # End-to-end processing pipeline
│   ├── 0.Functions.R            # QC helper functions (MAD cutoffs, Scrublet)
│   ├── 0.Preprocess.R           # Per-sample QC, gene filtering, doublet removal
│   ├── 1.Merge_datasets.R       # Merge per-sample Seurat objects → H5 / h5ad
│   ├── 2.Pegasus_Functions.py   # Pegasus QC helper module
│   ├── 2.Pegasus_run.py         # Pegasus clustering + QC
│   ├── 3.Manual_Major_Annotation.py  # Major cell-type annotation
│   ├── 4.Transfer_MAnnot_To_Seurat.R # Transfer Pegasus annotations to Seurat
│   ├── 5.Ref_SCANVI_Annot.py    # scVI/scANVI reference-based myeloid annotation
│   ├── 6.Metacells.R            # hdWGCNA metacell construction + Harmony
│   ├── 7.DEG_Metacells.R        # DEG: MTC_456 vs MTC_123
│   ├── 8.scDRS.py               # Disease relevance scoring (scDRS)
│   ├── 9.Liana.py               # Cell-cell interaction (LIANA)
│   ├── CCI_IREA_Association.py  # CCI × IREA enrichment association
│   └── utils.R                  # Shared color palettes, I/O helpers
│
├── Figures/
│   ├── Figure1/
│   │   ├── 1.Figure1_B_CCAmetadata.R
│   │   ├── 1.Figure1_C_UMAP.R
│   │   ├── 1.Figure1_D_CrumblR.R
│   │   └── 1.Figure1_E_Dreamlet.R
│   ├── Figure2/
│   │   ├── 2.Figure2_A_UMAP.R
│   │   ├── 2.Figure2_B_scRDS.py
│   │   ├── 2.Figure2_C_CrumblR.R
│   │   ├── 2.Figure2_D_CL_DEGs.R
│   │   └── 2.Figure2_E_GRN.R
│   ├── Figure3/
│   │   ├── 3.Figure3_A_Liana.py
│   │   ├── 3.Figure3_B_CellChat.R
│   │   └── 3.Figure3_C_IREA.R
│   ├── Figure4/
│   │   └── 4.Figure4_ExpTrends.R
│   ├── SFigure1/
│   │   └── SFigure1.R
│   ├── SFigure2/
│   │   ├── 1.Figure1D_CT_Markers.R
│   │   ├── SFigure2_A_Co_FreshMG_ICH.py
│   │   └── SFigure2_BCD.R
│   ├── SFigure3/
│   │   ├── SFigure_3A.R
│   │   ├── SFigure_3B.R
│   │   └── SFigure_3C.py
│   ├── SFigure4/
│   │   └── SFigure4.R
│   └── SFigure5/
│       └── SFigure5.py
│
├── Overview.jpeg                # Graphical overview of the study workflow
├── .zenodo.json                 # Zenodo metadata for code deposit
└── README.md                    # This file
```

---

## MTC groupings

| Group | Clusters | Outcome |
|-------|----------|---------|
| MTC_123 | MTC_1, MTC_2, MTC_3 | Good (mRS 0–3) |
| MTC_456 | MTC_4, MTC_5, MTC_6 | Poor (mRS 4–6) |
| MTC_Prolif | MTC_7 | Proliferating |

---

## Usage

All scripts accept command-line arguments — no hardcoded paths.

```bash
# 1. Per-sample QC
Rscript scripts/0.Preprocess.R \
    sample_name raw.h5 output_dir \
    TRUE TRUE 5 0.2 FALSE FALSE 15 FALSE \
    /path/to/MitoCarta.csv /path/to/gencode.gtf

# 2. Merge samples
Rscript scripts/1.Merge_datasets.R sample1.rds sample2.rds ... \
    metadata.csv output.rds output.h5Seurat output.h5ad

# 3. Pegasus QC + clustering
python scripts/2.Pegasus_run.py input.h5ad output_dir True

# 4. scANVI myeloid annotation
python scripts/5.Ref_SCANVI_Annot.py \
    reference.h5ad query.h5ad output.h5ad predictions.csv ./

# 5. Metacell construction
Rscript scripts/6.Metacells.R \
    input.rds Subclass FALSE 20 4 \
    output.pdf predictions.csv umap.csv output.rds

# 6. DEG: MTC_456 vs MTC_123
Rscript scripts/7.DEG_Metacells.R metacells.rds output_dir

# 7. Reproduce Figure 3C (IREA)
Rscript Figures/Figure3/3.Figure3_C_IREA.R \
    6a.IREA_MTC123.csv 6b.IREA_MTC456.csv Figure3C.pdf
```

---

## System requirements

**R packages (key)**

| Package | Version |
|---------|---------|
| Seurat | 4.9.9.9049 |
| dreamlet | 0.99.16 |
| variancePartition | 1.31.9 |
| crumblr | 0.99.6 |
| hdWGCNA | 0.2.19 |
| harmony | 0.1.1 |
| lme4 | ≥ 1.1 |
| lmerTest | ≥ 3.1 |
| emmeans | ≥ 1.8 |
| dittoSeq | 1.12.0 |
| clusterProfiler | 4.8.1 |
| SCopeLoomR | 0.13.0 |
| ComplexHeatmap | 2.16.0 |
| patchwork | 1.2.0 |
| archive | ≥ 1.1 |

**Python packages (key)**

| Package | Version |
|---------|---------|
| pegasuspy | ≥ 1.7 |
| scvi-tools | 1.0.4 |
| liana | 1.0.3 |
| scdrs | 1.0.3 |
| scanpy | 1.9.6 |
| anndata | 0.10.3 |
| scrublet | 0.2.3 |
| pandas | ≥ 2.0 |
| numpy | ≥ 1.23 |
| matplotlib | ≥ 3.7 |

---

## License

MIT License. See [LICENSE](LICENSE) for details.

Data files are released under CC BY 4.0 via the Zenodo data deposit (DOI: [DATA_DOI — to be added]).
