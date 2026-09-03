# Author: Dimitrios Kyriakis
import sys
import os
print(sys.version)
print(sys.executable)
import yaml
import sys
import warnings
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import scvi

import random

# Set a specific seed value
random.seed(94705)
scvi.settings.seed = 94705
TF_CPP_MIN_LOG_LEVEL=0

# ======================== READ DATA ==========================
print(sys.argv[0:])
ref_h5ad               = sys.argv[1]
query_h5ad             = sys.argv[2]
output_h5ad            = sys.argv[3]
output_csv_predictions = sys.argv[4]
path = sys.argv[5] if len(sys.argv) > 5 else 'result/5.FreshMG_Mapping/'

os.makedirs(path, exist_ok=True)
scvi_dir_path          = path + "FreshMG_ADAM_model_scvi_for_Biopsies_onlyMG/"
scanvi_dir_path_subtype  = path + "FreshMG_ADAM_model_scanvi_subtype_for_Biopsies_onlyMG/"
scanvi_dir_path_subclass = path + "FreshMG_ADAM_model_scanvi_subclass_for_Biopsies_onlyMG/"

# Load the source and target datasets
adata_ref   = sc.read(ref_h5ad)
adata_ref.obs['subclass']
adata_ref.layers['counts'] = adata_ref.layers['raw'].copy()

# Load the query datasets
adata_query = sc.read(query_h5ad)
adata_query.layers['counts'] = adata_query.X.copy()
adata_query.obs['Donor'] = adata_query.obs['donor']


# Select common genes
ref_genes = adata_ref.var_names
set2 = set(adata_query.var_names)
common_genes = [x for x in ref_genes if x in set2]
adata_ref = adata_ref[:, common_genes].copy()
adata_query = adata_query[:, common_genes].copy()


# Preprocess the datasets if needed
scvi.model.SCVI.setup_anndata(adata_ref, batch_key="Donor", layer="counts")
sc.pp.highly_variable_genes(adata_ref, n_top_genes=2000, batch_key="Donor", subset=True)
adata_query = adata_query[:, adata_ref.var_names].copy()


#%%  ================= Train the SCVI model on the source dataset =========================
if os.path.isdir(scvi_dir_path):
    print("Load SCVI REF model of Subtype")
    vae_ref = scvi.model.SCVI.load(scvi_dir_path, adata_ref)
else:
    scvi.model.SCVI.setup_anndata(adata_ref, batch_key="Donor", layer="counts")
    arches_params = dict(
        use_layer_norm="both",
        use_batch_norm="none",
        encode_covariates=True,
        dropout_rate=0.2,
        n_layers=2,
    )

    vae_ref = scvi.model.SCVI(adata_ref, **arches_params)
    vae_ref.train()
    vae_ref.save(scvi_dir_path, overwrite=True)
    adata_ref.obsm["X_scVI"] = vae_ref.get_latent_representation()

    adata_ref.obs["labels_scanvi"] = adata_ref.obs['subtype'].values
    vae_ref_scanvi = scvi.model.SCANVI.from_scvi_model(
        vae_ref,
        unlabeled_category="Unknown",
        labels_key="labels_scanvi",
    )
    vae_ref_scanvi.train(max_epochs=20, n_samples_per_label=100)
    adata_ref.obsm["X_scANVI"] = vae_ref_scanvi.get_latent_representation()
    sc.pp.neighbors(adata_ref, use_rep="X_scANVI")
    sc.tl.leiden(adata_ref)
    sc.tl.umap(adata_ref)
    sc.pl.umap(
        adata_ref,
        color=["subclass", "subtype"],
        frameon=False,
        ncols=1,
    )
    vae_ref_scanvi.save(scanvi_dir_path_subtype, overwrite=True)


#%%=======================================================================================
# ========================== Predict using SCANVI  subtype ===============================
scvi.model.SCANVI.prepare_query_anndata(adata_query, scanvi_dir_path_subtype)
vae_q = scvi.model.SCANVI.load_query_data(
    adata_query,
    scanvi_dir_path_subtype,
)
vae_q.train(
    max_epochs=50,
    plan_kwargs=dict(weight_decay=0.0),
    check_val_every_n_epoch=10,
)

adata_query.obsm["X_scANVI"] = vae_q.get_latent_representation()
adata_query.obs["predictions_subtype"] = vae_q.predict()


df = adata_query.obs.groupby(["CellType", "predictions_subtype"]).size().unstack(fill_value=0)
norm_df = df / df.sum(axis=0)
plt.figure(figsize=(8, 8))
_ = plt.pcolor(norm_df)
_ = plt.xticks(np.arange(0.5, len(df.columns), 1), df.columns, rotation=90)
_ = plt.yticks(np.arange(0.5, len(df.index), 1), df.index)
plt.xlabel("Predicted")
plt.ylabel("Observed")
plt.savefig(path+'MG_heatmap_predictions_CellType_vs_subtype.png')


#%%=======================================================================================
# ================= Train the SCANVI model on the source dataset (subclass) =========================
if os.path.isdir(scanvi_dir_path_subclass):
    print("Load SCANVI REF model of Subclass")
else:
    adata_ref.obs["labels_scanvi"] = adata_ref.obs['subclass'].values
    vae_ref_scanvi = scvi.model.SCANVI.from_scvi_model(
        vae_ref,
        unlabeled_category="Unknown",
        labels_key="labels_scanvi",
    )
    vae_ref_scanvi.train(max_epochs=20, n_samples_per_label=100)
    adata_ref.obsm["X_scANVI"] = vae_ref_scanvi.get_latent_representation()
    sc.pp.neighbors(adata_ref, use_rep="X_scANVI")
    sc.tl.leiden(adata_ref)
    sc.tl.umap(adata_ref)
    sc.pl.umap(
        adata_ref,
        color=["subclass", "subtype"],
        frameon=False,
        ncols=1,
    )
    vae_ref_scanvi.save(scanvi_dir_path_subclass, overwrite=True)


#%%=======================================================================================
# =========================== Predict using SCANVI subclass ==============================
scvi.model.SCANVI.prepare_query_anndata(adata_query, scanvi_dir_path_subclass)
vae_q = scvi.model.SCANVI.load_query_data(
    adata_query,
    scanvi_dir_path_subclass,
)
vae_q.train(
    max_epochs=50,
    plan_kwargs=dict(weight_decay=0.0),
    check_val_every_n_epoch=10,
)

adata_query.obsm["X_scANVI"] = vae_q.get_latent_representation()
adata_query.obs["predictions_subclass"] = vae_q.predict()
sc.pp.neighbors(adata_query, use_rep="X_scANVI")
sc.tl.leiden(adata_query)
sc.tl.umap(adata_query)


df = adata_query.obs.groupby(["CellType", "predictions_subclass"]).size().unstack(fill_value=0)
norm_df = df / df.sum(axis=0)
plt.figure(figsize=(8, 8))
_ = plt.pcolor(norm_df)
_ = plt.xticks(np.arange(0.5, len(df.columns), 1), df.columns, rotation=90)
_ = plt.yticks(np.arange(0.5, len(df.index), 1), df.index)
plt.xlabel("Predicted")
plt.ylabel("Observed")
plt.savefig(path+'MG_heatmap_predictions_CellType_vs_subclass.png')


#%%=======================================================================================
# ================================== Write Predictions ===================================
adata_query.obs[['predictions_subtype', 'predictions_subclass']].to_csv(output_csv_predictions, index=True)
adata_query.write_h5ad(output_h5ad)

output_csv_predictionsall = path+'MG_predictions_allmetadata.csv'
adata_query.obs.to_csv(output_csv_predictionsall, index=True)

output_umap = path+'MG_predictions_UMAP.csv'
df = pd.DataFrame(adata_query.obsm['X_umap'])
df.set_index(adata_query.obs_names)
df.columns = ['UMAP_1', 'UMAP_2']
df['CellName'] = adata_query.obs_names
df.to_csv(output_umap, index=True)

output_umap = path+'MG_predictions_scANVI.csv'
df = pd.DataFrame(adata_query.obsm['X_scANVI'])
df.set_index(adata_query.obs_names)
df['CellName'] = adata_query.obs_names
df.to_csv(output_umap, index=True)
