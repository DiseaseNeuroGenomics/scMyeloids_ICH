# Author: Dimitrios Kyriakis
# sc
import random
random.seed(123456)

import os
import sys
import pegasus as pg; import scanpy as sc;from scipy import stats;import numpy as np;import pandas as pd
# plotting
import matplotlib.pyplot as plt;
import seaborn as sns;

from scripts.PegasusFunctions import *

# ======================== READ DATA ==========================
print(sys.argv[0:])
input_file      = sys.argv[1]
remove_doublet  = sys.argv[2]
pca_regressed_harmony = sys.argv[3]
umap_coords     = sys.argv[4]
prefix          = sys.argv[5] if len(sys.argv) > 5 else 'result/2.Pegasus_run/'

os.makedirs(prefix, exist_ok=True)

run_doublet = False
data = pg.read_input(input_file)
# ================================================

marker_genes=['nCount_RNA','nFeature_RNA']

pg.violin(data, attrs=marker_genes, groupby='donor', dpi=300, wspace=1.2, hspace=1.2,bottom=1.2);plt.savefig(prefix+"1.QC.jpg")


pass_run="Run1"
file_exists = os.path.exists(prefix+pass_run+'_pre-doublet.zarr.zip')
if(file_exists):
    data = pg.read_input(prefix+pass_run+'_pre-doublet.zarr.zip')
else:
    data.obs.n_genes = data.obs.nFeature_RNA
    data.obs.n_counts =data.obs.nCount_RNA
    n_genes_lower, n_genes_upper = qc_boundary(data.obs.nFeature_RNA)
    n_counts_lower, n_counts_upper = qc_boundary(data.obs.nCount_RNA)
    print('genes lower: %d upper: %d' % (n_genes_lower, n_genes_upper))
    print('UMIs lower: %d upper: %d' % (n_counts_lower, n_counts_upper))

    plt.figure(figsize=(4, 4))
    sns.histplot(np.log10(data.obs.nCount_RNA));plt.axvline(np.log10(n_counts_lower), color='red');plt.axvline(np.log10(n_counts_upper), color='red');plt.savefig(prefix+"nCount.jpg")
    sns.displot(np.log10(data.obs.nFeature_RNA));plt.axvline(np.log10(n_genes_lower), color='red');plt.axvline(np.log10(n_genes_upper), color='red');plt.savefig(prefix+"nGenes.jpg")
    NODG = np.array(np.sum(data.X>0, axis=1))[:,0]
    temp = NODG.argsort()[::-1]
    ranks = np.empty_like(temp)
    ranks[temp] = np.arange(len(NODG))
    plt.figure(figsize=(6, 6));ax = sns.scatterplot(x=ranks,y=NODG);ax.set(yscale="log");plt.axhline(n_genes_lower, color='red');plt.axhline(n_genes_upper, color='red');plt.savefig(prefix+"nORG.jpg")

    NOU = np.array(np.sum(data.X, axis=1))[:,0]
    temp = NOU.argsort()[::1]
    ranks = np.empty_like(temp)
    ranks[temp] = np.arange(len(NOU))
    pg.qc_metrics(data,
                min_genes=n_genes_lower, max_genes=n_genes_upper,
                min_umis=n_counts_lower, max_umis=n_counts_upper,
                mito_prefix='MT-', percent_mito=10)
    pg.filter_data(data)
    plt.figure(figsize=(6, 6))
    ax = sns.scatterplot(x=ranks, y=NOU);ax.set(yscale="log");plt.axhline(n_counts_lower, color='red');plt.axhline(n_counts_upper, color='red');plt.savefig(prefix+"nODG.jpg")
    plt.figure(figsize=(4, 4))
    sns.scatterplot(x=data.obs.nFeature_RNA, y=data.obs.nCount_RNA, alpha=0.5, s=0.1)
    plt.axvline(n_genes_lower, color='red')
    plt.axvline(n_genes_upper, color='red')
    plt.axhline(n_counts_lower, color='red')
    plt.axhline(n_counts_upper, color='red');plt.savefig(prefix+"Boxed_Cells.jpg")
    pg.identify_robust_genes(data, percent_cells=0.05)
    data._inplace_subset_var(data.var['robust'])
    data

    data.var['non_ribosomal'] = [not (x.startswith("RPS") or x.startswith("RPL")) for x in data.var.index]
    data.var['non_mitochondrial'] = [not x.startswith("MT-") for x in data.var.index]
    data.var['robust'] = data.var['robust'] & data.var['non_ribosomal'] & data.var['non_mitochondrial']
    data.var[data.var['robust']]

    pg.log_norm(data)
    pg.write_output(data, prefix+'Cortex_norm.zarr.zip')

    pass_run="Run1"
    data = pass_run_fun(data, pass_run="Run1", prefix=prefix)

    pg.leiden(data, rep='pca_regressed_harmony', resolution=1.5)
    pg.write_output(data, prefix+pass_run+'_pre-doublet.zarr.zip')
    pg.scatter(data, attrs=['predicted_phase','percent.mt','donor','leiden_labels'], basis='umap', legend_loc='on data', dpi=150, wspace=0.1)
    plt.savefig(prefix+pass_run+"_Filtered_umap.png")

pg.scatter(data, attrs=['predicted_phase','percent.mt','donor','leiden_labels'], basis='umap', legend_loc='on data', dpi=150, wspace=0.1)
plt.savefig(prefix+pass_run+"_Filtered_umap_ledeng1.png")

pg.scatter(data, attrs=['is_doublet', 'scrublet_score'], basis='umap', legend_loc='on data', dpi=150, wspace=0.1)
plt.savefig(prefix+pass_run+"_Filtered_umap_scrublet.png")

pg.violin(data,
          attrs=['is_doublet', 'scrublet_score'],
          groupby='leiden_labels',
          dpi=300, wspace=1.2,
          hspace=1.2,bottom=1.2);
plt.savefig(prefix+"Vln_scrublet.jpg")

data.obs.is_doublet.value_counts()
data.obs['pred_dbl'] = [True if element==1 else False for element in data.obs.is_doublet ]
data.obs['pred_dbl'].value_counts()
data.obs['pred_dbl'][ data.obs.leiden_labels=='19'] = True
data.obs['pred_dbl'].value_counts()

pg.mark_doublets(data)
data.obs['demux_type'].value_counts()

pg.scatter(data, attrs='demux_type', basis='umap', dpi=150);plt.savefig(prefix+"Doublet_umap.png")


if (remove_doublet):
    pass_run="Run2"
    file_exists = os.path.exists(prefix+pass_run+'_post-doublet.zarr.zip')
    if(file_exists):
        data = pg.read_input(prefix+pass_run+'_post-doublet.zarr.zip')
    else:
        data.obs['demux_type'].value_counts()
        pg.qc_metrics(data, select_singlets=True)
        data.obs.passed_qc.value_counts()
        pg.filter_data(data)

        data = pass_run_fun(data, pass_run=pass_run, prefix=prefix)
        pg.write_output(data, prefix+pass_run+'_post-doublet.zarr.zip')
        pg.scatter(data, attrs=['predicted_phase','sex',
                                'donor','pred_dbl','leiden_labels'],
                            basis='umap', legend_loc='on data',
                            dpi=150, wspace=0.1);
        plt.savefig(prefix+pass_run+"_Metadata_umap.png")


pg.scatter(data, attrs=['predicted_phase','sex',
                        'donor','pred_dbl','leiden_labels'],
                        basis='umap', legend_loc='right margin',
                        dpi=150, wspace=0.1);
plt.savefig(prefix+pass_run+"_Metadata_umap_right_margin.png")

df = data.var['robust']
df.to_csv(prefix+'robust_genes.tsv',sep='\t')
# ===== PCA
data.obsm['X_pca_regressed_harmony']
df = pd.DataFrame(data.obsm['X_pca_regressed_harmony'])
list1 = list(range(1,df.shape[1]+1))
colnames_pc = list(map(lambda ls: "PC_" + str(ls), list1))
df.columns = colnames_pc
df.index = data.obs['orig.ident'].index
# ===== UMAP
df.to_csv(pca_regressed_harmony,sep='\t')
df = pd.DataFrame(data.obsm['X_umap'], columns = ['Umap_1','Umap_2'])
df.index = data.obs['orig.ident'].index
df.to_csv(umap_coords,sep='\t')
df.index = data.obs['orig.ident'].index

adata=data.to_anndata()
adata.write_h5ad(prefix+pass_run+'_post-doublet.h5ad')
