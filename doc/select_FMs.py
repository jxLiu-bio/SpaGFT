import SpaGFT.gft as gft
import SpaGFT.plot as plot
import numpy as np
import pandas as pd
import scanpy as sc
import os 
import time
import numpy as np
import random
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, ListedColormap

data_dir = "./datasets/151673"
adata = sc.read_visium(data_dir)
adata.raw = adata
sc.pp.filter_genes(adata, min_cells=10)
sc.pp.normalize_total(adata)
sc.pp.log1p(adata)

# Obtain SVGs under various ratio_low_freqs
n_top_genes = 1000
low_list = [0.5, 1, 2, 3, 4, 5, 8, 10]
svg_df = pd.DataFrame(index=low_list, columns=['svgs'])
for low in low_list: 
    gene_df = gft.rank_gene_smooth(adata,
                                   ratio_low_freq=low,
                                   ratio_high_freq=3,
                                   ratio_neighbors=1,
                                   spatial_info='spatial',
                                   filter_peaks=True,
                                   normalize_lap=False)
    gene_list = gene_df.index[:n_top_genes].tolist()
    svg_df.loc[low, 'svgs'] = gene_list

# Calculate overlap
low_index = ['ratio_' + str(i) for i in low_list]
overlap_df = pd.DataFrame(0, index=low_index, columns=low_index)
for ii in range(svg_df.shape[0]):
    for jj in range(ii, svg_df.shape[0]):
        overlap = np.intersect1d(svg_df.iloc[ii, 0], svg_df.iloc[jj, 0]).size
        overlap_df.iloc[ii, jj] = overlap
        overlap_df.iloc[jj, ii] = overlap

# Plot
plt.figure(figsize=(14, 12), dpi=300)
color_list = ['#FFFFFF', '#ec6666']
newcmap = LinearSegmentedColormap.from_list('chaos',color_list,N=256)
base_plot = sns.heatmap(overlap_df, annot=True, 
                        fmt="d", cmap=newcmap, vmin=850,
                        annot_kws={'size':20, 'weight':'bold','color':'white'})
plt.show()
plt.savefig("./results/sup_figs/heatmap_selection_fms_1000.png")

# freq_domain
freq_domain = adata.uns['fms_low']
freq_domain = pd.DataFrame(freq_domain, columns=['FM_' + str(i) for i in range(1, freq_domain.shape[1]+1)])
freq_domain[['x', 'y']] = adata.obsm['spatial']
freq_domain.to_csv("./results/sup_figs/151673_low_FMs.csv")