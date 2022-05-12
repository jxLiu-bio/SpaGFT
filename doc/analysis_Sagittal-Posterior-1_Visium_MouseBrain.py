import numpy as np
import pandas as pd
import scanpy as sc
import SpaGFT as spg
import os
import matplotlib.pyplot as plt
from plotnine import *


# Define data path
data_dir = "./datasets/"
sample = "Sagittal-Posterior-1_Visium_MouseBrain"
adata = sc.read_h5ad(data_dir + sample + "_expression.h5ad")
coord_df = pd.read_csv(data_dir + sample + "_xy.csv",
                       index_col=0)
adata.obs[['array_row', 'array_col']] = coord_df[['x', 'y']]
adata.var.index = adata.var._index
adata.raw.var.index = adata.raw.var._index
adata.var_names_make_unique()
sc.pp.filter_genes(adata, min_cells=10)
sc.pp.normalize_total(adata)
sc.pp.log1p(adata)


# *************************************
# Identify spatially variable gene
# *************************************
gene_df = spg.rank_gene_smooth(adata, 
                               spatial_info=['array_row', 'array_col'],
                               ratio_low_freq=0.5,
                               ratio_high_freq=3,
                               ratio_neighbors=1,
                               filter_peaks=True)
svg_list = gene_df[gene_df.cutoff_gft_score][gene_df.qvalue < 0.05].index.tolist()
gene_df = gene_df.loc[svg_list, :]
## genes in svg_list are spatially variable genes identified by SpaGFT

# Plot frequency signals of some genes and corresponding gene expressions
gene_signal = svg_list[:5]
for gene in gene_signal:
    # show the frequency signal
    freq_signal = adata[:, gene].varm['freq_domain_svg']
    freq_signal = np.ravel(freq_signal)
    plt.figure(figsize=(6, 2))
    low = list(range(adata.uns['fms_low'].shape[1]))
    plt.bar(low, freq_signal[low], color='#ca1c1c')
    high = list(range(len(low), freq_signal.size))
    plt.bar(high, freq_signal[high], color='#345591')
    ax = plt.gca()
    ax.set_ylabel("siganl")
    ax.spines['right'].set_color("none")
    ax.spines['top'].set_color("none")
    plt.ylim(0, 0.2)
    plt.xlim(0, freq_signal.size)
    plt.title("Gene: " + gene)
    plt.show()
    
    # show gene expression
    spg.scatter_gene_distri(adata, gene, shape='h',
                            spatial_info=['array_row', 'array_col'],
                            coord_ratio=1.3,
                            size=3)

# Visiualize SVG UMAP
spg.calculate_frequcncy_domain(adata,
                               ratio_low_freq=2,
                               ratio_high_freq=0,
                               ratio_neighbors=1,
                               spatial_info=['array_row', 'array_col'],
                               filter_peaks=False,
                               normalize_lap=False)
spg.plot.umap_svg(adata, size=2, coord_ratio=1.2)  


# ******************************
# Find tissue module
# ******************************
spg.gft.find_tissue_module(adata,
                           ratio_fms=2)
gene_df = adata.var.loc[svg_list, :]
# visualize TM
spg.plot.scatter_tm_expression(adata,
                               tm='0',
                               cmap='Spectral_r',
                               spatial_info=['array_row', 'array_col'],
                               coord_ratio=1.3,
                               size=3)
spg.plot.scatter_tm_expression(adata,
                               tm='0-0',
                               cmap='Spectral_r',
                               spatial_info=['array_row', 'array_col'],
                               coord_ratio=1.3,
                               size=3)
spg.plot.scatter_tm_binary(adata,
                            tm='0',
                            spatial_info=['array_row', 'array_col'],
                            coord_ratio=1.3,
                            size=3)

# Visiualize FMs
## low FMs
for i in range(1, 3):
    spg.plot.visualize_fms(adata, rank=i, size=3, coord_ratio=1.3)
## high FMs
for i in range(1, 3):
    spg.plot.visualize_fms(adata, rank=i, low=False, size=3, coord_ratio=1.3)

# results
results = adata.var.copy()
results.sort_values(by='gft_score', ascending=False, inplace=True)
results_folder ="./results/"
if not os.path.exists(results_folder):
    os.mkdir(results_folder)
results.to_csv(results_folder + sample + "_analysis.csv")

# **************************
# Enhancement
# **************************
new_adata = adata.copy()
spg.low_pass_enhancement(new_adata,
                        ratio_low_freq=15,
                        inplace=True)
for gene in svg_list[:5]:
    spg.scatter_gene_distri(new_adata, gene, shape='h',
                            spatial_info=['array_row', 'array_col'],
                            coord_ratio=1.3,
                            size=3)










