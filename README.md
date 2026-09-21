# Single-cell profiling of living human brain identifies myeloid states associated with six-month functional outcome after intracerebral hemorrhage

Analysis code, figures and supplementary tables for the manuscript submitted to **Nature Medicine**.

![Study overview: CD45+ cells were isolated from brain biopsies of 30 patients with acute intracerebral hemorrhage and profiled by single-cell RNA-seq. Six-month outcome was scored on the modified Rankin Scale and dichotomised into favorable (mRS 0-3) and unfavorable (mRS 4-6). The analysis arm comprises reference-based annotation, metacell aggregation, differential expression, gene regulatory network inference, cell-cell interaction analysis and computational drug repurposing, with prioritised compounds tested in a mouse ICH model.](Overview.jpeg)

<p align="center"><em>Study design and analysis workflow (manuscript Fig. 1a).</em></p>

---

## Authors

Dimitrios Kyriakis <sup>1,2,3,4,5 †</sup>,
Angelos Pavlopoulos <sup>6 †</sup>,
Xinyi Wang <sup>1,2,3,4 †</sup>,
Sarah Murphy <sup>1,2,3,4</sup>,
James M. Vicari <sup>1,2,3,4</sup>,
Rukmani Pandey <sup>7</sup>,
Gabriel E. Hoffman <sup>1,2,3,4,9</sup>,
Steve P. Kleopoulos <sup>1,2,3,4</sup>,
Stathis Argyriou <sup>1,2,3,4</sup>,
Zhiping Shao <sup>1,2,3,4</sup>,
Joon Ho Seo <sup>2,3</sup>,
Alexander Skupin <sup>5</sup>,
Nikolaos K. Robakis <sup>3,8</sup>,
Georgios Voloudakis <sup>1,2,3,4</sup>,
Anastasios Georgakopoulos <sup>3</sup>,
John F. Fullard <sup>1,2,3,4</sup>,
Donghoon Lee <sup>1,2,3,4</sup>,
Christopher P. Kellner,
Panos Roussos <sup>1,2,3,4,9</sup>

<sup>†</sup> These authors contributed equally to this work.

**Corresponding authors**
1. Panos Roussos, <panagiotis.roussos@mssm.edu>
2. Dimitrios Kyriakis, <dimitrios.kyriakis@mssm.edu>

### Affiliations

1. Friedman Brain Institute, Icahn School of Medicine at Mount Sinai, New York, NY, USA
2. Center for Disease Neurogenomics, Icahn School of Medicine at Mount Sinai, New York, NY, USA
3. Department of Psychiatry, Icahn School of Medicine at Mount Sinai
4. Department of Genetics and Genomic Sciences
5. Luxembourg Center for Systems Biomedicine (LCSB), University of Luxembourg, Esch-sur-Alzette, Luxembourg
6. Department of Pharmacology, Medical School of Athens, National and Kapodistrian University of Athens, Athens, Greece
7. Institute for Translational Medicine and Pharmacology, Icahn School of Medicine at Mount Sinai, New York, NY 10029, USA
8. Department of Neuroscience, Icahn School of Medicine at Mount Sinai
9. Mental Illness Research, Education and Clinical Centers, James J. Peters VA Medical Center, Bronx, New York

---

## Abstract

Intracerebral hemorrhage (ICH) causes high mortality and disability, but the human immune programs associated with recovery remain poorly defined. We profiled 93,378 CD45-positive cells from surgical brain biopsies obtained from 30 patients with acute ICH and examined cellular states to six-month functional outcome. Myeloid cells showed the largest outcome-associated transcriptional differences. Favorable outcome was associated with CX3CR1-rich surveillance and complement-related programs, whereas unfavorable outcome was associated with an SPP1-rich lipid-scavenging state with differential SPI1, PPARγ, and mTORC1-related activity. SPP1–integrin/CD44 signaling converged with findings from an independent human ICH cohort, and the outcome-associated myeloid axis mapped to stroke-associated programs in mice. Transcriptome-guided drug repurposing prioritized mTOR inhibition and retinoid signaling. In exploratory mouse studies, everolimus- and all-trans retinoic acid-treated groups showed better rotarod performance and altered spatial distributions of IBA1-positive myeloid morphologies. These findings define acute human myeloid states linked to recovery and nominate pathways for further preclinical evaluation.

---

## Repository structure

```
.
├── code/                          analysis source, mirrors the working tree
│   ├── 00_setup.R                 paths, constants, shared libraries
│   ├── 00_helpers_stats.R         statistical helper functions
│   ├── 00_theme_colors.R          ggplot theme and palettes
│   ├── build_supplementary_tables.R
│   ├── main/                      preprocessing and core pipeline
│   └── figures/
│       ├── Figure1/  Figure2/  Figure3/  Figure4/
│       ├── Supplementary/         Extended Data figures
│       └── qc/                    cross-check scripts
├── figures/                       final panels, PNG and PDF
│   ├── Figure_1..4
│   └── extended_data/Extended_Data_Figure_1..5
├── supplementary_tables/          Supplementary Tables 1-20, CSV
└── data/                          pointers to hosted data, see data/README.md
```

## Pipeline order

Run from `code/`. Scripts source `00_setup.R`, which defines `DERIVED_DIR`, `TABLE_DIR` and `FIG_DIR`.

| Stage | Scripts |
|---|---|
| Preprocessing and QC | `main/0.Preprocess.R`, `main/1.Merge_datasets.R`, `main/2.Pegasus_run.py` |
| Cell-type annotation | `main/3.Manual_Major_Annotation.py`, `main/4.Transfer_MAnnot_To_Seurat.R`, `main/5.Ref_SCANVI_Annot.py` |
| Metacell construction | `main/6.Metacells.R` |
| Disease relevance (scDRS) | `main/8.scDRS.py`, `main/9.scDRS_supplementary_table.py`, `figures/Supplementary/2026_scDRS.py` |
| Gene regulatory networks | `main/scenic_GRN_ExtendedGenes.lsf.sh`, `main/extract_loom_regulons.py`, `figures/Figure2/Figure2_G_*` |
| Cell-cell communication | `main/9.Liana.py`, `figures/Figure3/*` |
| Drug repurposing | see [voloudakislab/antagonist](https://github.com/voloudakislab/antagonist) |
| Mouse behaviour and morphology | `figures/Figure4/Figure4_Final_corrected_nums.R` |

### Which script makes which panel

| Panel | Script |
|---|---|
| Fig. 1b | `figures/Figure1/Figure1_B_CCA_heatmap.R` |
| Fig. 1c | `figures/Figure1/Figure1_C_umap.R` |
| Fig. 1d | `figures/Figure1/Figure1_D_mRS_effect_forest.R` |
| Fig. 1e | `figures/Figure1/Figure1_E_deg_summary_bar.R` |
| Fig. 1f | `figures/Figure1/Figure1_F_myeloid_enrichment.R` |
| Fig. 2a | `figures/Figure2/Figure2_A_umap_metacells.R` |
| Fig. 2b | `figures/Figure2/2.Figure2_B_scRDS.ipynb` |
| Fig. 2c | `figures/Figure2/Figure2_C_mRS_effect_tree.R` |
| Fig. 2d | `figures/Figure2/Figure2_D_heatmap_top_genes.R` |
| Fig. 2e | `figures/Figure2/Figure2_E_compareCluster_enrichment.R` |
| Fig. 2f | `figures/Figure2/Figure2_F_SAMC_composition_bar.R` |
| Fig. 2g | `figures/Figure2/Figure2_G_FINAL_dotplot.R` |
| Fig. 3a | `figures/Figure3/Figure3_A_irea_cytokine_plot.R` |
| Fig. 3b | `figures/Figure3/figure3b.py` |
| Fig. 3c | `figures/Figure3/figure3c.R` |
| Fig. 4b | `figures/Figure2/Figure2_H_hallmark_gsea_mTORC1.R` |
| Fig. 4c,d,e | `figures/Figure4/Figure4_Final_corrected_nums.R` |
| Extended Data Fig. 1 | `figures/Supplementary/SFigure1_donor_annotations.R` |
| Extended Data Fig. 2c,d,e,f | `figures/Supplementary/SFigure2_myeloid_subtypes.R` (batch wrapper: `SFigure2_myeloid_subtypes.lsf`). The same script writes Supplementary Table 6 and both source-data files |
| Extended Data Fig. 3b | `figures/Figure2/2.Figure2_B_scRDS.ipynb` |
| Extended Data Fig. 4 | `figures/Figure2/Figure2_G_02_network_degree.R` |
| Extended Data Fig. 5 | `figures/Figure4/Figure4_Final_corrected_nums.R` |

Image segmentation, morphological feature extraction and morphotype classification live in the separate repository [koniplus/microglia-zone-morphology](https://github.com/koniplus/microglia-zone-morphology).

---

## Supplementary tables

| # | Contents | Panel | File |
|---|---|---|---|
| 1 | Donor metadata | Extended Data Fig. 1 | `Supplementary_Table_01_donor_metadata_ExtDataFig1.csv` |
| 2 | CCA metadata | Fig. 1b | `Supplementary_Table_02_CCA_metadata_Fig1b.csv` |
| 3 | crumblr cell-type composition | Fig. 1d | `Supplementary_Table_03_crumblr_celltype_composition_Fig1d.csv` |
| 4 | DEGs associated with mRS, all immune cell types | Fig. 1e | `Supplementary_Table_04_DEGs_mRS_all_celltypes_Fig1e.csv` |
| 5 | Myeloid GSEA, full GO Biological Process | Fig. 1f | `Supplementary_Table_05_myeloid_GSEA_GOBP_Fig1f.csv` |
| 6 | dreamlet DEGs across myeloid subtypes | Extended Data Fig. 2e | `Supplementary_Table_06_dreamlet_DEGs_myeloid_subtypes_ExtDataFig2e.csv` |
| 7 | scDRS disease-relevance analysis | Fig. 2b | `Supplementary_Table_07_scDRS_disease_relevance_Fig2b.csv` |
| 8 | crumblr mRS associations across MTCs | Fig. 2c | `Supplementary_Table_08_crumblr_mRS_MTC_Fig2c.csv` |
| 9 | DEGs between MTC_456 and MTC_123 | Fig. 2d | `Supplementary_Table_09_DEGs_MTC456_vs_MTC123_Fig2d.csv` |
| 10 | GSEA, GO Molecular Function | Fig. 2e | `Supplementary_Table_10_GSEA_GOMF_Fig2e.csv` |
| 11 | Predicted cell-type composition | Fig. 2f | `Supplementary_Table_11_predicted_celltype_composition_Fig2f.csv` |
| 12 | Fisher TF enrichment analysis | Fig. 2g | `Supplementary_Table_12_fisher_TF_enrichment_Fig2g.csv` |
| 13 | Gene regulatory network degree ranking | Extended Data Fig. 4 | `Supplementary_Table_13_GRN_degree_ranking_ExtDataFig4.csv` |
| 14 | IREA results | Fig. 3a | `Supplementary_Table_14_IREA_results_Fig3a.csv` |
| 15 | LIANA cell–cell interaction results | Fig. 3b | `Supplementary_Table_15_LIANA_interactions_Fig3b.csv` |
| 16 | CellChat cell–cell interaction results | Fig. 3c | `Supplementary_Table_16_CellChat_interactions_Fig3c.csv` |
| 17 | Computational drug-repurposing compound results | — | `Supplementary_Table_17_drug_repurposing_compounds.csv` |
| 18 | Hallmark GSEA results | Fig. 4b | `Supplementary_Table_18_hallmark_GSEA_Fig4b.csv` |
| 19 | Behavioral analysis | Fig. 4c | `Supplementary_Table_19_behavioral_analysis_Fig4c.csv` |
| 20 | Morphology analysis | Fig. 4d,e | `Supplementary_Table_20_morphology_analysis_Fig4de.csv` |

### Source data

Not numbered supplementary tables. These are the full outputs behind two Extended Data panels, provided for transparency.

| Contents | Panel | File |
|---|---|---|
| crumblr myeloid-subtype composition: effect size, nominal p, and Benjamini–Hochberg FDR across the 13 subtypes | Extended Data Fig. 2d | `source_data/ExtDataFig2d_crumblr_subtype_composition.csv` |
| GO Molecular Function over-representation, complete results per subtype and direction | Extended Data Fig. 2f | `source_data/ExtDataFig2f_GOMF_enrichment.csv` |

---

## Software

Analyses ran in R and Python. Core packages, with the versions used:

**R** — Seurat, hdWGCNA v0.2.19, crumblr v0.99.6, dreamlet, variancePartition v1.31.9, muscat, Harmony v0.1.1, igraph, ggraph, clusterProfiler v4.8.3, org.Hs.eg.db v3.18.0, GO.db v3.18.0, metafor, lme4, emmeans

**Python** — Pegasus v1.7.0, scanpy v1.9.6, scDRS v1.0.3, pySCENIC v0.12.1, arboreto v0.1.6, LIANA v0.1.12, Cellpose v3.1.1, scikit-image, skan

`environment.yml` at the repository root recreates the conda environment used for the R analyses:

```bash
conda env create -f environment.yml -n ICH_STROKE_Env
```

Each script prints `sessionInfo()` on exit; the version block is written to the run log.

---

## Data availability

The single-cell dataset, clinical metadata and analysis outputs are available through
<https://cellxgene.cziscience.com/collections/d3c3e028-f91c-481e-9560-922fe94da67b>

See [`data/README.md`](data/README.md) for what is hosted where.

## Ethics approval

The human component of this study was approved by the Institutional Review Board of the Icahn School of Medicine at Mount Sinai (protocol number STUDY-18-01012A). Written informed consent was obtained from participants or their legally authorized representatives. All animal experiments were conducted according to institutional guidelines and were approved by the Institutional Animal Care and Use Committee of the Icahn School of Medicine at Mount Sinai (IACUC-2014-0271).

## Code availability

- Single-cell analysis: [DiseaseNeuroGenomics/scMyeloids_ICH](https://github.com/DiseaseNeuroGenomics/scMyeloids_ICH)
- Image-analysis pipeline: [koniplus/microglia-zone-morphology](https://github.com/koniplus/microglia-zone-morphology)
- Computational drug repurposing: [voloudakislab/antagonist](https://github.com/voloudakislab/antagonist), antagonist v1.0.0

## Competing interests

The authors declare no competing interests.

## Funding

This work was supported by the following: U.S. National Institutes of Health (NIH) NIA grants 2R56AG008200 (Nikolaos K. Robakis), 2R01AG008200 (Nikolaos K. Robakis), and 2R01NS047229 (Anastasios Georgakopoulos and Nikolaos K. Robakis); and the AP Slaner Family Award to Nikolaos K. Robakis. Additionally, this study was supported by the NIH, Bethesda, MD under award numbers R01AG067025, R01AG082185, R01AG065582, R01MH125246, R01AG050986, R01AG095776, U24AG087563, RF1NS147363, and U19AG097398. This work was supported in part through the computational and data resources and staff expertise provided by Scientific Computing and Data at the Icahn School of Medicine at Mount Sinai and supported by the Clinical and Translational Science Award (CTSA) grant UL1TR004419 from the National Center for Advancing Translational Sciences.

## Author contributions

DK, AS, DL, CK and PR conceived and designed the study. CK performed the endoscopic hematoma-evacuation surgeries and provided the brain biopsy specimens together with the clinical, radiographic and six-month outcome data. SK, SA, ZS and JF performed the tissue dissociation, FACS enrichment and single-cell RNA-sequencing library preparation. DK, AS, DL and PR designed the single-cell analysis strategy. DK performed the single-cell analysis, including quality control and preprocessing, cell-type and myeloid-subtype annotation, compositional and differential-expression analyses, metacell construction and clustering, gene-regulatory-network inference, disease-relevance scoring, cytokine-response and cell–cell communication analyses, and all data visualization. GV performed the computational drug repurposing. AP, JV and AG performed the mouse ICH experiments. AP and AG performed the behavioral assessments. JS contributed to the design of the mouse experimental setup. RP prepared the mouse brain samples for immunostaining. XW and SM performed the immunofluorescence staining and image acquisition. XW performed the image segmentation, morphological feature extraction and unsupervised morphotype classification. DK and GH performed the compositional and statistical analysis of the morphotype and behavioral data. DK, AP, DL and PR wrote the manuscript. DL and PR supervised the study. All authors contributed to the article and approved the submitted version.

## Citation

Kyriakis D et al. Single-cell profiling of living human brain identifies myeloid states associated with six-month functional outcome after intracerebral hemorrhage. *Submitted*, 2026.

## License

Code is released under the MIT License, see [LICENSE](LICENSE). Figures, supplementary tables and derived data accompany the manuscript and are subject to the publisher's terms.
