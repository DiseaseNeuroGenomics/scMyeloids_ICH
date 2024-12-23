import numpy as np
import pegasus as pg; import scanpy as sc;from scipy import stats;
import numpy as np;import pandas as pd
import matplotlib.pyplot as plt;

def qc_boundary(counts, k=3):
    x = np.log1p(counts)
    mad = stats.median_abs_deviation(x)
    return np.exp(np.median(x)-k*mad), np.exp(np.median(x)+k*mad)



def calc_MG_subtype_scores(pg_data):
    pg.calc_signature_score(pg_data, {"MG_homeostatic": ['P2RY12', 'TMEM119', 'GPR34', 'JUN',
                                                    'OLFML3', 'CSF1R', 'HEXB', 'MERTK',
                                                    'RHOB', 'CX3CR1', 'TGFBR1', 'TGFB1',
                                                    'MEF2A', 'MAFB', 'SALL1', 'EGR1',
                                                    'SOCS3','OLFML3','SPI1', 'SMAD3', 'ENTPD1']})
    pg.calc_signature_score(pg_data, {"MG_inflammatory": ['SPP1', 'ITGAX', 'AXL', 'LILRB4', 'CLEC7A', 'CCL2', 'CSF1']})
    pg.calc_signature_score(pg_data, {"MG_chemokine": ['CCL2', 'CCL3','CCL4']})
    pg.calc_signature_score(pg_data, {"MG_coreTF": ['MEF2A', 'MAFB', 'SALL1', 'EGR1']})
    pg.calc_signature_score(pg_data, {"M0": ['RUNX1','SPI1','IRF8']})
    pg.calc_signature_score(pg_data, {"MG_cytoskeleton": ['TMSB4X', 'ACTB', 'TMSB10', 'EPB41L2', 'MACF1', 'SRGAP2', 'UTRN']})
    pg.calc_signature_score(pg_data, {"MG_amyloid": ['ITM2B', 'CST3']})
    pg.calc_signature_score(pg_data, {"MG_repair": ['ATM', 'ERCC6']})
    pg.calc_signature_score(pg_data, {"MG_membrane": ['SORL1', 'SLCO2B1']})
    pg.calc_signature_score(pg_data, {"MG_DAM": ['AXL', 'APOE', 'CLEC7A', 'ITGAX', 'LGALS3', 'CST7']}) # https://www.nature.com/articles/s41583-018-0057-5
    pg.calc_signature_score(pg_data, {"MGnD": ['BHLHE40', 'TFLEC', 'ETS2', 'ATF3']})
    return(pg_data)



def pass_run_fun(pg_data,pass_run):
    prefix='result/Cortex/2.Pegasus_run/'
    adata=pg_data.to_anndata()
    # find highly variable genes
    sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5, batch_key='donor')
    print(adata.var.highly_variable.value_counts())
    sc.pl.highly_variable_genes(adata);plt.savefig(prefix+pass_run+"_sc_hvg.jpg")
    # remove non-robust genes from hvg
    adata.var.loc[adata.var.robust==False, 'highly_variable'] = False
    print(adata.var.highly_variable.value_counts())
    # copy
    pg_data.var.highly_variable_features = adata.var.highly_variable
    print(pg_data.var.highly_variable_features.value_counts())
    ### cell cycle gene score based on [Tirosh et al. 2015 | https://science.sciencemag.org/content/352/6282/189]
    pg.calc_signature_score(pg_data, 'cell_cycle_human')
    pg_data.obs['predicted_phase'].value_counts()
    pg_data.obs['CC_diff'] = pg_data.obs['G1/S'] - pg_data.obs['G2/M']
    ### pca/harmony/umap
    pg.pca(pg_data, n_components=15)
    pg.regress_out(pg_data, attrs=['nCount_RNA','percent.mt','CC_diff'])
    pg.run_harmony(pg_data, batch='donor', rep='pca_regressed', max_iter_harmony=20)
    # plt.savefig(prefix+"_haarmony_converged.jpg")
    print("Calculate neighbors")
    pg.neighbors(pg_data, rep='pca_regressed_harmony', use_cache=False)#, dist='cosine')
    pg.umap(pg_data, rep='pca_regressed_harmony') # rep='pca_regressed_harmony'
    ### clustering
    pg.leiden(pg_data, rep='pca_regressed_harmony', resolution=1.5)
    return(pg_data)