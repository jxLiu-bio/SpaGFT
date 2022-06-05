import pandas as pd
from plotnine import *
import scipy.sparse as ss
import numpy as np
import SpaGFT.gft as gft
import matplotlib.pyplot as plt
import scanpy as sc

def svg_freq_signal(adata, gene,
                    domain='freq_domain_svg',
                    figsize=(6, 2),
                    dpi=100, 
                    colors=['#CA1C1C', '#345591'],
                    return_fig=False, **kwargs):
    # Show the frequency signal
    freq_signal = adata[:, gene].varm[domain]
    freq_signal = np.ravel(freq_signal)
    plt.figure(figsize=figsize, dpi=dpi)
    low = list(range(adata.uns['fms_low'].shape[1]))
    plt.bar(low, freq_signal[low], color=colors[0])
    high = list(range(len(low), freq_signal.size))
    plt.bar(high, freq_signal[high], color=colors[1])
    ax = plt.gca()
    ax.set_ylabel("siganl")
    ax.spines['right'].set_color("none")
    ax.spines['top'].set_color("none")
    plt.ylim(0, 0.2)
    plt.xlim(0, freq_signal.size)
    plt.title("Gene: " + gene)
    plt.show()
    if return_fig:
        return ax

def tm_freq_signal(adata, tm,
                  domain='freq_domain_svg', figsize=(6, 2),
                  dpi=100, color='#CA1C1C',
                  y_range=[0, 0.08],
                  return_fig=False, **kwargs):
    # Show the frequency signal
    freq_signal = signal = adata.uns['freq_signal_tm'].loc[tm, :].values
    plt.figure(figsize=figsize, dpi=dpi)
    low = list(range(len(freq_signal)))
    plt.bar(low, freq_signal, color=color)
    plt.grid(False)
    ax = plt.gca()
    ax.set_ylabel("siganl")
    ax.spines['right'].set_color("none")
    ax.spines['top'].set_color("none")
    plt.ylim(y_range[0], y_range[1])
    plt.xlim(0, freq_signal.size)
    plt.title(tm)
    plt.show()
    if return_fig:
        return ax

def subTm_freq_signal(adata, subTM,
                  domain='freq_domain_svg', figsize=(6, 2),
                  dpi=100, color='#CA1C1C',
                  y_range=[0, 0.08],
                  return_fig=False, **kwargs):
    # Show the frequency signal
    freq_signal = signal = adata.uns['freq_signal_subTM'].loc[subTM, :].values
    plt.figure(figsize=figsize, dpi=dpi)
    low = list(range(len(freq_signal)))
    plt.bar(low, freq_signal, color=color)
    ax = plt.gca()
    plt.grid(False)
    ax.set_ylabel("siganl")
    ax.spines['right'].set_color("none")
    ax.spines['top'].set_color("none")
    plt.ylim(y_range[0], y_range[1])
    plt.xlim(0, freq_signal.size)
    plt.title(subTM)
    plt.show()
    if return_fig:
        return ax
    
def gene_signal_umap(adata, svg_list, colors=['#C81E1E', '#BEBEBE'], 
                     return_fig=False, n_neighbors=15, 
                     random_state=0, **kwargs):
    freq_domain = adata.varm['freq_domain_svg']
    freq_domain = pd.DataFrame(freq_domain)
    freq_domain.index = adata.var_names
    umap_adata = sc.AnnData(freq_domain)
    sc.pp.neighbors(umap_adata, n_neighbors=n_neighbors, use_rep='X')
    sc.tl.umap(umap_adata, random_state=0)
    adata.varm['freq_umap_svg'] = umap_adata.obsm['X_umap']
    print("""The umap coordinates of genes when identify SVGs could be found in 
          adata.varm['freq_umap_svg']""")
    # svg_list
    umap_adata.obs['SpaGFT'] = 'Non-SVGs'
    umap_adata.obs.loc[svg_list, 'SpaGFT'] = 'SVGs'
    umap_adata.obs['SpaGFT'] = pd.Categorical(umap_adata.obs['SpaGFT'], 
                                              categories=['SVGs', 'Non-SVGs'],
                                              ordered=True)
    umap_adata.uns['SpaGFT_colors'] = colors
    t = sc.pl.umap(umap_adata, color="SpaGFT", return_fig=return_fig, 
                   **kwargs)



def scatter_gene_distri(adata, 
                        gene, 
                        size=3, 
                        shape='h',
                        cmap='magma',
                        spatial_info=['in_tissue', 'array_row'],
                        coord_ratio=0.7, return_plot=False):
    if ss.issparse(adata.X):
        plot_df = pd.DataFrame(adata.X.todense(), index=adata.obs_names,
                           columns=adata.var_names)
    else:
        plot_df = pd.DataFrame(adata.X, index=adata.obs_names,
                           columns=adata.var_names)
    if spatial_info in adata.obsm_keys():
        plot_df['x'] = adata.obsm[spatial_info][:, 0]
        plot_df['y'] = adata.obsm[spatial_info][:, 1]
    elif set(spatial_info) <= set(adata.obs.columns):
        plot_coor = adata.obs
        plot_df = plot_df[gene]
        plot_df = pd.DataFrame(plot_df)
        plot_df['x'] = plot_coor.loc[:, spatial_info[0]].values
        plot_df['y'] = plot_coor.loc[:, spatial_info[1]].values
    plot_df['radius'] = size
    plot_df = plot_df.sort_values(by=gene, ascending=True)
    for i in [gene]:
        base_plot = (ggplot() + geom_point(plot_df, aes(x='x', y='y', fill=i), 
                                           shape=shape, stroke=0.1, size=size) +
                      xlim(min(plot_df.x)-1, max(plot_df.x)+1) + 
                      ylim(min(plot_df.y)-1, max(plot_df.y)+1) + 
                      scale_fill_cmap(cmap_name=cmap) + 
                      coord_equal(ratio=coord_ratio) +
                      theme_classic() +
                      theme(legend_position=('right'),
                            legend_background=element_blank(),
                            legend_key_width=4,
                            legend_key_height=50)
                      )
        print(base_plot)
    if return_plot:
        return base_plot

def umap_svg_cluster(adata,
                     svg_list=None,
                     size=3,
                     shape='o',
                     cmap='magma', 
                     spatial_info=['array_row', 'array_col'],
                     coord_ratio=1, 
                     return_plot=True):
    

    if svg_list == None:
        tmp_df = adata.var.copy()
        svg_list = tmp_df[tmp_df.cutoff_gft_score][tmp_df.qvalue < 0.05].index
    plot_df = adata.varm['gft_umap_tm']
    plot_df = pd.DataFrame(plot_df)
    plot_df.index = adata.var_names
    plot_df.columns = ['UMAP_1', 'UMAP_2']

    plot_df.loc[svg_list, 'gene'] = 'SVG'
    plot_df['gene'] = pd.Categorical(plot_df['gene'], 
                                     categories=['SVG', 'Non-SVG'],
                                     ordered=True)
    plot_df['radius'] = size
    # plot
    base_plot = (ggplot(plot_df, aes(x='UMAP_1', y='UMAP_2', fill='gene')) + 
                 geom_point(size=size, color='white', stroke=0.25) +
                 scale_fill_manual(values=colors) +
                 theme_classic() +
                 coord_equal(ratio=coord_ratio))
    print(base_plot)
    if return_plot:
        return base_plot
    

def scatter_tm_expression(adata, tm, size=3, shape='o', cmap='magma',
                          spatial_info=['array_row', 'array_col'],
                          coord_ratio=0.7, return_plot=False):
    if '-' in tm:
        tm = 'tm-' + tm.split('-')[0] + "_subTm-" + tm.split('-')[1]
        plot_df = adata.obsm['subTm_expression']
    else:
        tm = 'tm_' + tm
        plot_df = adata.obsm['tm_expression']
    plot_df = pd.DataFrame(plot_df)
    if spatial_info in adata.obsm_keys():
        plot_df = plot_df[tm]
        plot_df['x'] = adata.obsm[spatial_info][:, 0]
        plot_df['y'] = adata.obsm[spatial_info][:, 1]
    elif set(spatial_info) <= set(adata.obs.columns):
        plot_coor = adata.obs
        plot_df = plot_df[tm]
        plot_df = pd.DataFrame(plot_df)
        plot_df['x'] = plot_coor.loc[:, spatial_info[0]].values
        plot_df['y'] = plot_coor.loc[:, spatial_info[1]].values
    plot_df['radius'] = size
    plot_df = plot_df.sort_values(by=tm, ascending=True)
    base_plot = (ggplot() + geom_point(plot_df, aes(x='x', y='y', fill=tm), 
                                       shape=shape, stroke=0.1, size=size) +
                  xlim(min(plot_df.x)-1, max(plot_df.x)+1) + 
                  ylim(min(plot_df.y)-1, max(plot_df.y)+1) + 
                  scale_fill_cmap(cmap_name=cmap) + 
                  coord_equal(ratio=coord_ratio) +
                  theme_classic() +
                  theme(legend_position=('right'),
                        legend_background=element_blank(),
                        legend_key_width=4,
                        legend_key_height=50)
                  )
    print(base_plot)
    if return_plot:
        return base_plot
    
def scatter_tm_binary(adata, tm, size=3, shape='h',
                          spatial_info=['array_row', 'array_col'],
                          colors=['#CA1C1C','#CCCCCC'],
                          coord_ratio=0.7, return_plot=False):
    if '-' in tm:
        tm = 'tm-' + tm.split('-')[0] + "_subTm-" + tm.split('-')[1]
        plot_df = adata.obsm['subTm_region']
    else:
        tm = 'tm_' + tm
        plot_df = adata.obsm['tm_region']
    plot_df = pd.DataFrame(plot_df)
    if spatial_info in adata.obsm_keys():
        plot_df['x'] = adata.obsm[spatial_info][:, 0]
        plot_df['y'] = adata.obsm[spatial_info][:, 1]
    elif set(spatial_info) <= set(adata.obs.columns):
        plot_coor = adata.obs
        plot_df = plot_df[tm]
        plot_df = pd.DataFrame(plot_df)
        plot_df['x'] = plot_coor.loc[:, spatial_info[0]].values
        plot_df['y'] = plot_coor.loc[:, spatial_info[1]].values
    plot_df['radius'] = size
    plot_df[tm] = plot_df[tm].values.astype(int)
    plot_df[tm] = plot_df[tm].values.astype(str)
    plot_df[tm] = pd.Categorical(plot_df[tm], 
                                   categories=['1', '0'],
                                   ordered=True)
    base_plot = (ggplot() + geom_point(plot_df, aes(x='x', y='y', fill=tm), 
                                       shape=shape, stroke=0.1, size=size) +
                  xlim(min(plot_df.x)-1, max(plot_df.x)+1) + 
                  ylim(min(plot_df.y)-1, max(plot_df.y)+1) + 
                  scale_fill_manual(values=colors) + 
                  coord_equal(ratio=coord_ratio) +
                  theme_classic() +
                  theme(legend_position=('right'),
                        legend_background=element_blank(),
                        legend_key_width=4,
                        legend_key_height=50)
                  )
    print(base_plot)
    if return_plot:
        return base_plot

def umap_svg(adata, svg_list=None, colors=['#CA1C1C','#CCCCCC'], size=2,
             coord_ratio=0.7, return_plot=False):
    if 'gft_umap_svg' not in adata.varm_keys():
        raise KeyError("Please run SpaGFT.calculate_frequcncy_domain firstly")
    plot_df = adata.varm['gft_umap_svg']
    plot_df = pd.DataFrame(plot_df)
    plot_df.index = adata.var_names
    plot_df.columns = ['UMAP_1', 'UMAP_2']
    plot_df['gene'] = 'Non-SVG'
    if svg_list == None:
        tmp_df = adata.var.copy()
        svg_list = tmp_df[tmp_df.cutoff_gft_score][tmp_df.qvalue < 0.05].index
    plot_df.loc[svg_list, 'gene'] = 'SVG'
    plot_df['gene'] = pd.Categorical(plot_df['gene'], 
                                     categories=['SVG', 'Non-SVG'],
                                     ordered=True)
    plot_df['radius'] = size
    # plot
    base_plot = (ggplot(plot_df, aes(x='UMAP_1', y='UMAP_2', fill='gene')) + 
                 geom_point(size=size, color='white', stroke=0.25) +
                 scale_fill_manual(values=colors) +
                 theme_classic() +
                 coord_equal(ratio=coord_ratio))
    print(base_plot)
    if return_plot:
        return base_plot
    
def visualize_fms(adata, rank=1, low=True, size=3, cmap='magma',
                 spatial_info=['array_row', 'array_col'], shape='h',
                 coord_ratio=0.7, return_plot=False):
    if low == True:
        plot_df = pd.DataFrame(adata.uns['fms_low'])
        plot_df.index = adata.obs.index
        plot_df.columns = ['low_FM_' + str(i + 1) for i in range(plot_df.shape[1])]
        if spatial_info in adata.obsm_keys():
            plot_df['x'] = adata.obsm[spatial_info][:, 0]
            plot_df['y'] = adata.obsm[spatial_info][:, 1]
        elif set(spatial_info) <= set(adata.obs.columns):
            plot_coor = adata.obs
            plot_df = plot_df['low_FM_' + str(rank)]
            plot_df = pd.DataFrame(plot_df)
            plot_df['x'] = plot_coor.loc[:, spatial_info[0]].values
            plot_df['y'] = plot_coor.loc[:, spatial_info[1]].values
        plot_df['radius'] = size
        base_plot = (ggplot() + geom_point(plot_df, aes(x='x', y='y',
                                                   fill='low_FM_' + str(rank)), 
                                           shape=shape, stroke=0.1, size=size) +
                      xlim(min(plot_df.x)-1, max(plot_df.x)+1) + 
                      ylim(min(plot_df.y)-1, max(plot_df.y)+1) + 
                      scale_fill_cmap(cmap_name=cmap) + 
                      coord_equal(ratio=coord_ratio) +
                      theme_classic() +
                      theme(legend_position=('right'),
                            legend_background=element_blank(),
                            legend_key_width=4,
                            legend_key_height=50)
                      )
        print(base_plot)
        
    else:
        plot_df = pd.DataFrame(adata.uns['fms_high'])
        plot_df.index = adata.obs.index
        plot_df.columns = ['high_FM_' + str(i + 1) for i in\
                           range(adata.uns['fms_high'].shape[1])]
        if spatial_info in adata.obsm_keys():
            plot_df['x'] = adata.obsm[spatial_info][:, 0]
            plot_df['y'] = adata.obsm[spatial_info][:, 1]
        elif set(spatial_info) <= set(adata.obs.columns):
            plot_coor = adata.obs
            plot_df = plot_df['high_FM_' + str(plot_df.shape[1] - rank + 1)]
            plot_df = pd.DataFrame(plot_df)
            plot_df['x'] = plot_coor.loc[:, spatial_info[0]].values
            plot_df['y'] = plot_coor.loc[:, spatial_info[1]].values
        plot_df['radius'] = size
        base_plot = (ggplot() + geom_point(plot_df, aes(x='x', y='y',
                            fill='high_FM_' + \
                                str(adata.uns['fms_high'].shape[1] - rank + 1)), 
                            shape=shape, stroke=0.1, size=size) +
                      xlim(min(plot_df.x)-1, max(plot_df.x)+1) + 
                      ylim(min(plot_df.y)-1, max(plot_df.y)+1) + 
                      scale_fill_cmap(cmap_name=cmap) + 
                      coord_equal(ratio=coord_ratio) +
                      theme_classic() +
                      theme(legend_position=('right'),
                            legend_background=element_blank(),
                            legend_key_width=4,
                            legend_key_height=50)
                      )
        print(base_plot)
        
    if return_plot:
        return base_plot
    
    
    
    
    
  