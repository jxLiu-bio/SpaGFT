import SpaGFT.gft as gft
import SpaGFT.plot as plot
import numpy as np
import pandas as pd
import scanpy as sc
import os 
import time
import numpy as np
import random


# Set data path
data_path = "D:/project/SpaGFT/datasets"
sample = "HE_Visum_MouseBrain"

# ***************************************
# Load data
# ***************************************
adata = sc.read_h5ad(os.path.join(data_path, sample + "_expression.h5ad"))
adata.raw = adata
coord = pd.read_csv(os.path.join(data_path, sample + "_xy.csv"),
                          index_col=0)
adata.var.index = adata.var._index.tolist()
adata.raw.var.index = adata.raw.var._index.tolist()
# coord = coord.reindex(adata.obs_names)
# adata.obs[coord.columns] = coord
adata.obs["x"] = coord.values[:, 0]
adata.obs["y"] = coord.values[:, 1]
adata.var_names_make_unique()
sc.pp.filter_genes(adata, min_cells=10)
sc.pp.normalize_total(adata)
sc.pp.log1p(adata)

# ****************************************
# Find svgs
# ****************************************
gene_df= gft.rank_gene_smooth(adata,
                             ratio_low_freq=0.5,
                             ratio_high_freq=3,
                             ratio_neighbors=1,
                             filter_peaks=True,
                             spatial_info=['x', 'y'],
                             normalize_lap=False,
                             S=5)
svg_list = gene_df[gene_df.cutoff_gft_score]\
    [gene_df.qvalue < 0.05].index.tolist()
    
print(len(svg_list))


# **************************************
# Obtain frequency domain
# **************************************
gft.calculate_frequcncy_domain(adata,
                               spatial_info=['x', 'y'], 
                               return_freq_domain=False)
gft.find_tissue_module(adata,
                       spatial_info=['x', 'y'],
                       ratio_fms=2,
                       ratio_neighbors=1,
                       resolution=2)

# **************************************
# imputation
# **************************************
gft.low_pass_imputation(adata,
                        spatial_info=['x', 'y'])


























# # ****************************************
# # Find tissue modules
# # ****************************************
# sc.settings.set_figure_params(dpi=200, frameon=False, figsize=(7, 7),
#                               facecolor='white')
# adata = raw.copy()
# sc.pp.log1p(adata)
# gft.calculate_frequcncy_domain(adata,
#                                    ratio_low_freq=4,
#                                    ratio_high_freq=0,
#                                    ratio_neighbors=1,
#                                    spatial_info=['x', 'y'],
#                                    filter_peaks=False,
#                                    normalization_lap=False)
# gft_adata = sc.AnnData(adata.varm['freq_domain'])
# gft_adata.obs_names = adata.var_names
# gft_adata.obs['spatially_variable_genes'] = 'Non-SVG'
# gft_adata.obs.loc[svg_list, 'spatially_variable_genes'] = 'SVG'
# gft_adata.obs['spatially_variable_genes'] = \
#     pd.Categorical(gft_adata.obs['spatially_variable_genes'],
#                                                            categories=['Non-SVG', 'SVG'],
#                                                            ordered=True)
# sc.pp.neighbors(gft_adata, use_rep='X')
# sc.tl.umap(gft_adata,)
# # sc.pl.umap(gft_adata, color='spatially_variable_genes', s=30)
# gft_adata.uns['spatially_variable_genes_colors'] = ["#C0C0C0", "#C81E1E"]
# sc.pl.umap(gft_adata, color='spatially_variable_genes', s=30)

# # clustering
# tmp_adata = gft_adata[svg_list, :]
# sc.pp.neighbors(tmp_adata, use_rep='X')
# sc.tl.louvain(tmp_adata, resolution=3)
# all_tms = np.unique(tmp_adata.obs.louvain)
# gft_adata.obs['tissue_module_genes'] = 'None'
# gft_adata.obs.loc[svg_list, 'tissue_module_genes'] = tmp_adata.obs.louvain
# sc.pl.umap(gft_adata, color='tissue_module_genes', s=30)
# gft_adata.uns['tissue_module_genes_colors'][-1] = "#C0C0C0"
# sc.pl.umap(gft_adata, color='tissue_module_genes', s=30)

# # for tm in all_tms:
# #     for gene in gft_adata.obs[gft_adata.obs.tissue_module_genes == tm].index[-10:]:
# #         plot.scatter_gene_distri(adata, gene=gene, ratio=0.7, size=3)
# #     for gene in gft_adata.obs[gft_adata.obs.tissue_module_genes == tm].index[:10]:
# #         plot.scatter_gene_distri(adata, gene=gene, ratio=0.7, size=3)
# #     sc.pl.umap(gft_adata, color='tissue_module_genes', s=30)
    
# results = pd.DataFrame(gft_adata.obsm['X_umap'], index=gft_adata.obs_names,
#                       columns=['UMAP_1', 'UMAP_2'])
# results['tissue_module'] = gft_adata.obs['tissue_module_genes'] 
# results.to_csv("D:\\project\\SpaGFT\\results\\" + sample + "low-4_noFilter.csv")
# freq_domain = pd.DataFrame(adata.varm['freq_domain'],
#                             index=adata.var_names)
# freq_domain.columns = ["FM_" + str(i+1) for i in range(freq_domain.shape[1])]
# freq_domain.to_csv("D:\\project\\SpaGFT\\results\\" + sample + "low-4_freqDomain.csv")



