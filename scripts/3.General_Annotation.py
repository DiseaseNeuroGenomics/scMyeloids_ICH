import random
random.seed(123456)
import os
import sys
import pegasus as pg; 
import scanpy as sc;
from scipy import stats;
import numpy as np;
import pandas as pd
# plotting
import matplotlib.pyplot as plt;
import seaborn as sns;


# data
prefix='ICH_Stroke/result/3.Annotation/'
from PegasusFunctions import *
run_doublet = False
os.chdir('ICH_Stroke/')
print(os.getcwd())


# ======================== READ DATA ==========================
print(sys.argv[1:])
input_file = sys.argv[1]
activation_score_genes = sys.argv[2]
output_zip = sys.argv[3]
output_metadata = sys.argv[4]
output_h5ad = sys.argv[5]


data = pg.read_input(input_file)

# =================================== Auto Annotation =================================
# =================================== Auto Annotation =================================
pg.de_analysis(data, cluster='leiden_labels');marker_dict = pg.markers(data);pg.write_results_to_excel(marker_dict,prefix+'df_leiden.xlsx')

# Annotate based on  markers
celltype_dict = pg.infer_cell_types(data, markers = 'ICH_project/Data/human_brain_immune_cell_markers.json')
cluster_names = pg.infer_cluster_names(celltype_dict)
pg.annotate(data, name='anno_dh', based_on='leiden_labels', anno_dict=cluster_names)

pg.scatter(data, attrs=['donor','leiden_labels','anno_dh'], basis='umap', legend_loc='on data', dpi=150, wspace=0.1)
plt.savefig(prefix+"Test.png")
pg.dendrogram(data, rep='pca_regressed_harmony', groupby='leiden_labels', panel_size=(5,5), dpi=100)
plt.savefig(prefix+"Dendrogram.png")
# ------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------

# ============================ Manual Annotation =================================
# ============================ Manual Annotation =================================
cluster_names = [
"Microglia",# 1
"Microglia",# 2
"Microglia",# 3
"Microglia",# 4
"Microglia",# 5
"Microglia",# 6
"Microglia",# 7
"Microglia",# 8
"Microglia",# 9

"Neutrophil",# 10
"Microglia",# 11
"Prolif_Mono",# 12
"T Cells",# 13
"Astrocytes",# 14
"Oligodendrocytes",# 15
"Mural",# 16
"Memory B cells"]# 17

# "Microglia",# 10
# "Microglia",# 13

# "Microglia",# 20
# "Microglia"]# 21

pg.annotate(data, name='class', based_on='leiden_labels', anno_dict=cluster_names)
pg.scatter(data, attrs=['leiden_labels','class','anno_dh'], basis='umap', legend_loc='on data', dpi=150, wspace=0.1)
plt.savefig(prefix+"Annotations.png")
pg.dendrogram(data, rep='pca_regressed_harmony', groupby='class', panel_size=(5,5), dpi=100)
plt.savefig(prefix+"Dendrogram_Class.png")
# # ------------------------------------------------------------------------------------
# # ------------------------------------------------------------------------------------



# # ======================================= Activation Score =============================================================
# # ======================================= Activation Score =============================================================
df = pd.read_csv( activation_score_genes,header=0,usecols=["Module", "GeneSynbol"])
df = df.dropna()
adata=data.to_anndata()
ASCORE = {}
for i in np.unique(df.Module):
    print(i)
    module_genes = list(map(lambda x: x.upper(), df.GeneSynbol[df.Module==i].tolist()))
    list1_as_set = set(module_genes)
    intersection = list1_as_set.intersection(adata.raw.var_names)
    intersection_as_list = list(intersection)
    ASCORE[i]=intersection_as_list

sc.pl.matrixplot(adata, ASCORE, 'leiden_labels', dendrogram=True, cmap='Blues', standard_scale='var', colorbar_title='column scaled\nexpression')
plt.savefig(prefix+"Enzymatic_Heat_subclass.jpg")
pg.calc_signature_score(data, ASCORE)
# # ------------------------------------------------------------------------------------
# # ------------------------------------------------------------------------------------

adata=data.to_anndata()
sc.pl.scatter(adata,x="Micro/Myeloid Shared Act. Score",y="Microglial Identity Score", legend_loc='on data',color="leiden_labels")
plt.savefig(prefix+"Scatter_activation_vs_identity.jpg")



# # ======================================= PLOTS =============================================================
# # ======================================= PLOTS =============================================================
marker_genes_dict={'Ast': ['CAMK2G','CKB'],
    'Pro-inflammatory': ['HSPB1','FOSB','JUN','PLK2'],
    'Homeostatic': ['CD74','EEF1B2','GPR34','ALOX5AP','APOC2','APOE','C1QB','C1QC','C3'],
    'Cytoskeleton': ['NPIPB5'],
    'Histone': ['HIST1H2AE','HIST1H2BD','HIST1H2BF','HIST1H2BG','HIST1H2BJ','HIST1H3H'],
    'Interferon': ['IFI44L','IFIT2','MX2'],
    'Anti-inflammatory': ['GPNMB','SELENOP'],
    'Cell-cycle': ['STMN1'],
    'Chemokine': ['SPP1','CCL3','CCL4','CCL4L2','IL1B','EGR1'],
    'Mono': ['LGALS1','LRRK2','LYZ','S100A6','S100A8','S100A9','VIM'],
    'Mural': ['DSTN','SELENOM','TPM1'],'Olig': ['CNP','KIF1B','SLC44A1','TMEM144'],
    'T-cell': ['CD48','CXCR4','ETS1','GNLY','LY6E','UTRN']}
adata=data.to_anndata()
sc.pl.matrixplot(adata, marker_genes_dict, 'leiden_labels', dendrogram=True, cmap='Blues', standard_scale='var', colorbar_title='column scaled\nexpression')
plt.savefig(prefix+"All_Heat_subclass.jpg")
sc.pl.matrixplot(adata, marker_genes_dict, 'class', dendrogram=True, cmap='Blues', standard_scale='var', colorbar_title='column scaled\nexpression')
plt.savefig(prefix+"All_Heat_class.jpg")


pg.compo_plot(data, 'mRS','class', panel_size=(4, 1))
plt.savefig(prefix+"mRS_class.jpg")


marker_genes_dict={'Ast': ['GJA1'],
    'Microglia': ["P2RY12"],
    'Olig': ["PLP1"],
    #'Endo': ["LY6C1"],
    'T-cell': ['CD3D'],
    'B-cell': ['CD79A'],
    'Cell-cycle': ['STMN1'],
    'Chemokine': ['CCL3','CCL4','IL1B'],
    'Mono': ['VIM'],
    'Mural': ['TPM1']}
    'Pro-inflammatory': ['JUN'],
    'Cell-cycle': ['STMN1'],
    'Chemokine': ['CCL3','CCL4','IL1B']}


plt.figure(figsize=(15, 6))
adata=data.to_anndata()
adata = adata[adata.obs['class'].isin(['Astrocytes', 'Monocytes', 'Oligodendrocytes', 'T Cells','Memory B cells'])]
sc.tl.dendrogram(adata,groupby="class")
ax = sc.pl.tracksplot(adata, marker_genes_dict, groupby='class', dendrogram=False,figsize =(15,6))
plt.savefig(prefix+"Bars_expression.jpg")

# # ------------------------------------------------------------------------------------
# # ------------------------------------------------------------------------------------

df = pd.DataFrame(data.obs)
df.to_csv(output_metadata,sep='\t')


pg.write_output(data,output_h5ad,file_type="h5ad")
pg.write_output(data,output_zip)
