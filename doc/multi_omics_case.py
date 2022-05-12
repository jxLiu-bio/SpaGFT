import SpaGFT.gft as gft
import SpaGFT.plot as plot
import numpy as np
import pandas as pd
import scanpy as sc
import os 
import time
import numpy as np
import random

# data_dir = "./datasets/multi_omics"
# sample = "0713_DBiT-Seq_MouseEmbryos"
# adata = sc.read(os.path.join(data_dir, sample + "_expression.h5ad"))
# adata.raw = adata
# adata.var.index = adata.var._index
# adata.raw.var.index = adata.raw.var._index
# coord = pd.read_csv(os.path.join(data_dir, sample + "_location.csv"))
# coord.index = adata.obs_names
# adata.obs = coord
data_dir = "/bmbl_data/yuzhou/spatial_data/merge_data/"
sample = "0713_DBiT-Seq_MouseEmbryos"
adata = sc.read(os.path.join(data_dir + "expression", sample + "_expression.h5ad"))
adata.raw = adata
adata.var.index = adata.var._index
adata.raw.var.index = adata.raw.var._index
coord = pd.read_csv(os.path.join(data_dir + "location", sample + "_location.csv"))
coord.index = adata.obs_names
adata.obs = coord

# Split RNAs proteins
feature_list = adata.var._index.tolist()
protein_list = []
for i in feature_list:
    if '-pro' in i:
        protein_list.append(i)
rna_list = np.setdiff1d(feature_list, protein_list)

rna_adata = adata[:, rna_list]
pro_adata = adata[:, protein_list]

# Preprocessing
sc.pp.filter_genes(rna_adata, min_cells=10)
sc.pp.normalize_total(rna_adata)
sc.pp.log1p(rna_adata)
sc.pp.log1p(pro_adata)
# Find significant svgs
gene_df = gft.rank_gene_smooth(rna_adata,
                               ratio_low_freq=0.5,
                               ratio_high_freq=3,
                               ratio_neighbors=1,
                               filter_peaks=True,
                               normalize_lap=False,
                               spatial_info = ['x', 'y'])
svg_list = gene_df[gene_df.cutoff_gft_score][gene_df.qvalue < 0.05].index.tolist()

# calculate freq_domain
adata = adata[:, np.union1d(protein_list, svg_list)]
sc.pp.log1p(adata)
freq_df = gft.calculate_frequcncy_domain(adata,
                                         ratio_low_freq=2,
                                         ratio_high_freq=0,
                                         ratio_neighbors=1,
                                         filter_peaks=False,
                                         spatial_info=['x', 'y'],
                                         normalize_lap=False)
# clustering
gft_adata = sc.AnnData(freq_df.transpose())
gft_adata.obs['attribute'] = "rna"
gft_adata.obs.loc[protein_list, 'attribute'] = 'protein'
sc.pp.neighbors(gft_adata, use_rep='X')
sc.tl.umap(gft_adata)
for res in [0.5, 1, 1.5, 2]:
    sc.tl.louvain(gft_adata, resolution=res, random_state=0)
    sc.pl.umap(gft_adata, color='louvain')
    cross_df = gft_adata.obs.copy()
    cross_df = cross_df.sort_values(by='attribute')
    cross_df.to_csv("/bmbl_data/yuzhou/Python_results/" + sample+ "_low-2_res-" + \
                    str(res) + ".csv")
freq_df = freq_df.loc[:, cross_df.index]
freq_df = freq_df.transpose()
freq_df.to_csv("/bmbl_data/yuzhou/Python_results/" + sample+ "_low-2_freqDomain.csv")





