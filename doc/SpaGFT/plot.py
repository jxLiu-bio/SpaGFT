# -*- coding: utf-8 -*-

import pandas as pd
from plotnine import *
import scipy.sparse as ss
import numpy as np

def scatter_gene_distri(adata, gene, size=5, shape='o', color='auto',
                        spatial_names=['x', 'y']):
    if ss.issparse(adata.X):
        plot_df = pd.DataFrame(adata.X.todense(), index=adata.obs_names,
                           columns=adata.var_names)
    else:
        plot_df = pd.DataFrame(adata.X, index=adata.obs_names,
                           columns=adata.var_names)
    plot_coor = adata.obs
    plot_df = plot_df[gene]
    plot_df = pd.DataFrame(plot_df)
    plot_df['x'] = plot_coor.loc[:, spatial_names[0]].values
    plot_df['y'] = plot_coor.loc[:, spatial_names[1]].values
    plot_df['radius'] = size
    plot_df = plot_df.sort_values(by=gene, ascending=True)
    for i in [gene]:
        if color == 'auto':
            base_plot = (ggplot() + geom_point(plot_df, aes(x='x', y='y', fill=i), 
                                               shape=shape, stroke=0.1, size=size) +
                          xlim(min(plot_df.x)-1, max(plot_df.x)+1) + 
                          ylim(min(plot_df.y)-1, max(plot_df.y)+1) + 
                          scale_fill_cmap(cmap_name='magma') + 
                          coord_equal(ratio=1) +
                          theme_classic() +
                          theme(legend_position=('right'),
                                legend_background=element_blank(),
                                legend_key_width=4,
                                legend_key_height=50)
                          )
        else:
            base_plot = (ggplot() + geom_point(plot_df, aes(x='x', y='y', size=i,fill=i), 
                                                shape=shape, stroke=0.1, color='black') +
                          xlim(min(plot_df.x)-1, max(plot_df.x)+1) + 
                          ylim(min(plot_df.y)-1, max(plot_df.y)+1) +
                          scale_fill_manual(color) +
                          scale_size_area(min_size=0.3) +
                          theme_classic()
                          )
            
            
        print(base_plot)
    
    return base_plot
    
def scatter_group(plot_array, data_inf, axis_label=['Dim_1', 'Dim_2'],
                  data_label = 'group', shape='o', color=None, axis_range=None):
    plot_array = pd.DataFrame(plot_array, columns=axis_label)
    plot_array[data_label] = data_inf
    plot_array[data_label] = pd.Categorical(plot_array[data_label])
    if color == None:
        base_plot = (ggplot(plot_array, aes(x=axis_label[0], y=axis_label[1],
                            fill=data_label))
                     + geom_point()
                     + scale_fill_hue(s=0.99, l=0.65, h=0.0417, 
                                      color_space='husl'))
    else:
        base_plot = (ggplot(plot_array, aes(x=axis_label[0], y=axis_label[1],
                            fill=data_label))
                     + geom_point()
                     + scale_fill_manual(color))
        if not axis_range == None:
            base_plot = (ggplot(plot_array, aes(x=axis_label[0], y=axis_label[1],
                                fill=data_label))
                         + geom_point()
                         + scale_fill_manual(color)
                         + scale_x_continuous(limits=(axis_range[0],
                                                           axis_range[1]))
                         + scale_y_continuous(limits=(axis_range[2],
                                                           axis_range[3])))
    print(base_plot)

def pca_spatial_domain(adata, n_top=1000):
    import scanpy as sc
    tmp_adata = adata.copy().T
    tmp_adata = tmp_adata[tmp_adata.obs.sort_values(by='SVG_Rank').index, :]
    if not 'spatially_variable' in tmp_adata.obs_names:
        tmp_adata.obs['spatially_variable'] = 0
        tmp_adata.obs.iloc[:n_top, len(tmp_adata.obs.columns)-1] = 1
        tmp_adata.obs['spatially_variable'] = \
            pd.Categorical(tmp_adata.obs['spatially_variable'], 
                           categories=[1, 0])
    tmp_adata = tmp_adata[tmp_adata.obs.sort_values(by='SVG_Rank', 
                                                    ascending=False).index, :]
    sc.pp.pca(tmp_adata, n_comps=50)
    
    # plot
    scatter_group(plot_array=tmp_adata.obsm['X_pca'][:, :2],
                  data_inf=tmp_adata.obs['spatially_variable'].tolist(),
                  axis_label=['PC_1', 'PC_2'],
                  data_label='SVG',
                  color=['#00AFBB', '#FC4E07'])
    scatter_group(plot_array=tmp_adata[:n_top, :].obsm['X_pca'][:, :2],
                  data_inf=tmp_adata[:n_top, :].obs['spatially_variable'].tolist(),
                  axis_label=['PC_1', 'PC_2'],
                  data_label='SVG',
                  color=['#FC4E07'])
    scatter_group(plot_array=tmp_adata[n_top:, :].obsm['X_pca'][:, :2],
                  data_inf=tmp_adata[n_top:, :].obs['spatially_variable'].tolist(),
                  axis_label=['PC_1', 'PC_2'],
                  data_label='SVG',
                  color=['#00AFBB'])

def umap_spatial_domain(adata, n_top=1000):
    import scanpy as sc
    tmp_adata = adata.copy().T
    tmp_adata = tmp_adata[tmp_adata.obs.sort_values(by='SVG_Rank').index, :]
    if not 'spatially_variable' in tmp_adata.obs_names:
        tmp_adata.obs['spatially_variable'] = "Non"
        tmp_adata.obs.iloc[:n_top, len(tmp_adata.obs.columns)-1] = 1
    tmp_adata.obs.loc[gene_cluster_info.index, 'spatially_variable'] = gene_cluster_info.louvain
    tmp_adata = tmp_adata[tmp_adata.obs.sort_values(by='SVG_Rank', 
                                                    ascending=False).index, :]
    sc.pp.pca(tmp_adata, n_comps=50)
    obs = tmp_adata.obs
    tmp_adata = sc.AnnData(tmp_adata.obsm['X_pca'])
    tmp_adata.obs = obs
    sc.pp.neighbors(tmp_adata, use_rep='X')
    sc.tl.umap(tmp_adata)
    
    # plot
    scatter_group(plot_array=tmp_adata.obsm['X_umap'][:, :2],
                  data_inf=tmp_adata.obs['spatially_variable'].tolist(),
                  axis_label=['UMAP_1', 'UMAP_2'],
                  data_label='SVG',
                  color=['#00AFBB', '#FC4E07'])
    scatter_group(plot_array=tmp_adata[:n_top, :].obsm['X_umap'][:, :2],
                  data_inf=tmp_adata[:n_top, :].obs['spatially_variable'].tolist(),
                  axis_label=['UMAP_1', 'UMAP_2'],
                  data_label='SVG',
                  color=['#FC4E07'])
    scatter_group(plot_array=tmp_adata[n_top:, :].obsm['X_umap'][:, :2],
                  data_inf=tmp_adata[n_top:, :].obs['spatially_variable'].tolist(),
                  axis_label=['UMAP_1', 'UMAP_2'],
                  data_label='SVG',
                  color=['#00AFBB'])

def tsne_spatial_domain(adata, n_top=1000):
    import scanpy as sc
    tmp_adata = adata.copy().T
    tmp_adata = tmp_adata[tmp_adata.obs.sort_values(by='SVG_Rank').index, :]
    if not 'spatially_variable' in tmp_adata.obs_names:
        tmp_adata.obs['spatially_variable'] = 0
        tmp_adata.obs.iloc[:n_top, len(tmp_adata.obs.columns)-1] = 1
        tmp_adata.obs['spatially_variable'] = \
            pd.Categorical(tmp_adata.obs['spatially_variable'], 
                           categories=[1, 0])
    tmp_adata = tmp_adata[tmp_adata.obs.sort_values(by='SVG_Rank', 
                                                    ascending=False).index, :]
    sc.pp.pca(tmp_adata, n_comps=50)
    obs = tmp_adata.obs
    tmp_adata = sc.AnnData(tmp_adata.obsm['X_pca'])
    tmp_adata.obs = obs
    sc.pp.neighbors(tmp_adata)
    sc.tl.tsne(tmp_adata)
    
    # plot
    scatter_group(plot_array=tmp_adata.obsm['X_tsne'][:, :2],
                  data_inf=tmp_adata.obs['spatially_variable'].tolist(),
                  axis_label=['tSNE_1', 'tSNE_2'],
                  data_label='SVG',
                  color=['#00AFBB', '#FC4E07'])
    scatter_group(plot_array=tmp_adata[:n_top, :].obsm['X_tsne'][:, :2],
                  data_inf=tmp_adata[:n_top, :].obs['spatially_variable'].tolist(),
                  axis_label=['tSNE_1', 'tSNE_2'],
                  data_label='SVG',
                  color=['#FC4E07'])
    scatter_group(plot_array=tmp_adata[n_top:, :].obsm['X_tsne'][:, :2],
                  data_inf=tmp_adata[n_top:, :].obs['spatially_variable'].tolist(),
                  axis_label=['tSNE_1', 'tSNE_2'],
                  data_label='SVG',
                  color=['#00AFBB'])

def fm_frequency_domain(adata, n_top=1000):
    import scanpy as sc
    tmp_adata = sc.AnnData(adata.varm['freq_domain'])
    tmp_adata.obs = adata.var
    tmp_adata = tmp_adata[tmp_adata.obs.sort_values(by='SVG_Rank').index, :]
    if not 'spatially_variable' in tmp_adata.obs_names:
        tmp_adata.obs['spatially_variable'] = "Non"
        tmp_adata.obs.iloc[:n_top, len(tmp_adata.obs.columns)-1] = 1
    tmp_adata.obs.loc[gene_cluster_info.index, 'spatially_variable'] = gene_cluster_info.louvain
    
    tmp_adata = tmp_adata[tmp_adata.obs.sort_values(by='SVG_Rank', 
                                                    ascending=False).index, :]
    
    scatter_group(plot_array=tmp_adata.X[:, 1:3],
                  data_inf=tmp_adata.obs['spatially_variable'].tolist(),
                  axis_label=['FM_1', 'FM_2'],
                  data_label='SVG',
                  color=['#00AFBB', '#FC4E07'])
    scatter_group(plot_array=tmp_adata[:n_top, :].X[:, 1:3],
                  data_inf=tmp_adata[:n_top, :].obs['spatially_variable'].tolist(),
                  axis_label=['FM_1', 'FM_2'],
                  data_label='SVG',
                  color=['#FC4E07'])
    scatter_group(plot_array=tmp_adata[n_top:, :].X[:, 1:3],
                  data_inf=tmp_adata[n_top:, :].obs['spatially_variable'].tolist(),
                  axis_label=['FM_1', 'FM_2'],
                  data_label='SVG',
                  color=['#00AFBB'])

def fm_umap_spectral_domain(adata, n_top=1000):
    import scanpy as sc
    tmp_adata = sc.AnnData(adata.varm['freq_domain'])
    tmp_adata.obs = adata.var
    tmp_adata = tmp_adata[tmp_adata.obs.sort_values(by='SVG_Rank').index, :]
    if not 'spatially_variable' in tmp_adata.obs_names:
        tmp_adata.obs['spatially_variable'] = "Non"
        tmp_adata.obs.iloc[:n_top, len(tmp_adata.obs.columns)-1] = 1
    tmp_adata.obs.loc[gene_cluster_info.index, 'spatially_variable'] = gene_cluster_info.louvain
    tmp_adata = tmp_adata[tmp_adata.obs.sort_values(by='SVG_Rank', 
                                                    ascending=False).index, :]
    
    sc.pp.pca(tmp_adata, n_comps=50)
    sc.pp.neighbors(tmp_adata)
    sc.tl.umap(tmp_adata)
    
    scatter_group(plot_array=tmp_adata.obsm['X_umap'][:, :2],
                  data_inf=tmp_adata.obs['spatially_variable'].tolist(),
                  axis_label=['UMAP_1', 'UMAP_2'],
                  data_label='SVG',
                  color=['#00AFBB', '#FC4E07'])
    scatter_group(plot_array=tmp_adata[:n_top, :].obsm['X_umap'][:, :2],
                  data_inf=tmp_adata[:n_top, :].obs['spatially_variable'].tolist(),
                  axis_label=['UMAP_1', 'UMAP_2'],
                  data_label='SVG',
                  color=['#FC4E07'])
    scatter_group(plot_array=tmp_adata[n_top:, :].obsm['X_umap'][:, :2],
                  data_inf=tmp_adata[n_top:, :].obs['spatially_variable'].tolist(),
                  axis_label=['UMAP_1', 'UMAP_2'],
                  data_label='SVG',
                  color=['#00AFBB'])

def fm_tsne_spectral_domain(adata, n_top=1000):
    import scanpy as sc
    tmp_adata = sc.AnnData(adata.varm['freq_domain'])
    tmp_adata.obs = adata.var
    tmp_adata = tmp_adata[tmp_adata.obs.sort_values(by='SVG_Rank').index, :]
    if not 'spatially_variable' in tmp_adata.obs_names:
        tmp_adata.obs['spatially_variable'] = 0
        tmp_adata.obs.iloc[:n_top, len(tmp_adata.obs.columns)-1] = 1
        tmp_adata.obs['spatially_variable'] = \
            pd.Categorical(tmp_adata.obs['spatially_variable'], 
                           categories=[0, 1])
    tmp_adata = tmp_adata[tmp_adata.obs.sort_values(by='SVG_Rank', 
                                                    ascending=False).index, :]
    
    sc.pp.pca(tmp_adata, n_comps=50)
    sc.pp.neighbors(tmp_adata)
    sc.tl.tsne(tmp_adata)
    
    scatter_group(plot_array=tmp_adata.obsm['X_tsne'][:, :2],
                  data_inf=tmp_adata.obs['spatially_variable'].tolist(),
                  axis_label=['tSNE_1', 'tSNE_2'],
                  data_label='SVG',
                  color=['#00AFBB', '#FC4E07'])
    scatter_group(plot_array=tmp_adata[:n_top, :].obsm['X_tsne'][:, :2],
                  data_inf=tmp_adata[:n_top, :].obs['spatially_variable'].tolist(),
                  axis_label=['tSNE_1', 'tSNE_2'],
                  data_label='SVG',
                  color=['#FC4E07'])
    scatter_group(plot_array=tmp_adata[n_top:, :].obsm['X_tsne'][:, :2],
                  data_inf=tmp_adata[n_top:, :].obs['spatially_variable'].tolist(),
                  axis_label=['tSNE_1', 'tSNE_2'],
                  data_label='SVG',
                  color=['#00AFBB'])

def fm_umap_direct(adata, n_top=1000):
    import scanpy as sc
    tmp_adata = sc.AnnData(adata.varm['freq_domain'][:, :50])
    tmp_adata.obs = adata.var
    tmp_adata = tmp_adata[tmp_adata.obs.sort_values(by='SVG_Rank').index, :]
    if not 'spatially_variable' in tmp_adata.obs_names:
        tmp_adata.obs['spatially_variable'] = 0
        tmp_adata.obs.iloc[:n_top, len(tmp_adata.obs.columns)-1] = 1
        tmp_adata.obs['spatially_variable'] = \
            pd.Categorical(tmp_adata.obs['spatially_variable'], 
                           categories=[0, 1])
    tmp_adata = tmp_adata[tmp_adata.obs.sort_values(by='SVG_Rank', 
                                                    ascending=False).index, :]
    sc.pp.neighbors(tmp_adata, use_rep='X')
    sc.tl.umap(tmp_adata)
    
    scatter_group(plot_array=tmp_adata.obsm['X_umap'][:, :2],
                  data_inf=tmp_adata.obs['spatially_variable'].tolist(),
                  axis_label=['UMAP_1', 'UMAP_2'],
                  data_label='SVG',
                  color=['#00AFBB', '#FC4E07'])
    scatter_group(plot_array=tmp_adata[:n_top, :].obsm['X_umap'][:, :2],
                  data_inf=tmp_adata[:n_top, :].obs['spatially_variable'].tolist(),
                  axis_label=['UMAP_1', 'UMAP_2'],
                  data_label='SVG',
                  color=['#FC4E07'],
                  axis_range=[-12, 3, -4, 4])
    scatter_group(plot_array=tmp_adata[n_top:, :].obsm['X_umap'][:, :2],
                  data_inf=tmp_adata[n_top:, :].obs['spatially_variable'].tolist(),
                  axis_label=['UMAP_1', 'UMAP_2'],
                  data_label='SVG',
                  color=['#00AFBB'])
    
def fm_tsne_direct(adata, n_top=1000):
    import scanpy as sc
    tmp_adata = sc.AnnData(adata.varm['freq_domain'][:, :50])
    tmp_adata.obs = adata.var
    tmp_adata = tmp_adata[tmp_adata.obs.sort_values(by='SVG_Rank').index, :]
    if not 'spatially_variable' in tmp_adata.obs_names:
        tmp_adata.obs['spatially_variable'] = 0
        tmp_adata.obs.iloc[:n_top, len(tmp_adata.obs.columns)-1] = 1
        tmp_adata.obs['spatially_variable'] = \
            pd.Categorical(tmp_adata.obs['spatially_variable'], 
                           categories=[0, 1])
    tmp_adata = tmp_adata[tmp_adata.obs.sort_values(by='SVG_Rank', 
                                                    ascending=False).index, :]
    sc.pp.neighbors(tmp_adata, use_rep='X')
    sc.tl.tsne(tmp_adata)
    
    scatter_group(plot_array=tmp_adata.obsm['X_tsne'][:, :2],
                  data_inf=tmp_adata.obs['spatially_variable'].tolist(),
                  axis_label=['tSNE_1', 'tSNE_2'],
                  data_label='SVG',
                  color=['#00AFBB', '#FC4E07'])
    scatter_group(plot_array=tmp_adata[:n_top, :].obsm['X_tsne'][:, :2],
                  data_inf=tmp_adata[:n_top, :].obs['spatially_variable'].tolist(),
                  axis_label=['tSNE_1', 'tSNE_2'],
                  data_label='SVG',
                  color=['#FC4E07'])
    scatter_group(plot_array=tmp_adata[n_top:, :].obsm['X_tsne'][:, :2],
                  data_inf=tmp_adata[n_top:, :].obs['spatially_variable'].tolist(),
                  axis_label=['tSNE_1', 'tSNE_2'],
                  data_label='SVG',
                  color=['#00AFBB'])
    


    
    
    
            
    
    

    

    
    

