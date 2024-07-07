import pandas as pd
from plotnine import *
import scipy.sparse as ss
import numpy as np
import scanpy as sc
from sklearn import preprocessing
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.gridspec import GridSpec
import math
import os
from matplotlib.colors import ListedColormap


def gene_freq_signal(adata,
                     gene,
                     domain='freq_domain_svg',
                     figsize=(6, 2),
                     dpi=100,
                     colors=['#CA1C1C', '#345591'],
                     save_path=None,
                     return_fig=False,
                     **kwargs):
    if isinstance(gene, str):
        freq_signal = adata[:, gene].varm[domain]
        freq_signal = np.ravel(freq_signal)
        plt.figure(figsize=figsize, dpi=dpi)
        low = list(range(adata.uns['identify_svg_data']['fms_low'].shape[1]))
        plt.bar(low, freq_signal[low], color=colors[0])
        high = list(range(len(low), freq_signal.size))
        plt.bar(high, freq_signal[high], color=colors[1])
        ax = plt.gca()
        ax.set_ylabel("signal")
        ax.spines['right'].set_color("none")
        ax.spines['top'].set_color("none")
        y_max = max(freq_signal)
        plt.ylim(0, y_max * 1.1)
        plt.xlim(0, freq_signal.size)
        plt.title("Gene: " + gene)
        plt.show()
        if save_path is not None:
            plt.savefig(f"{save_path}")
        if return_fig:
            return ax
    elif isinstance(gene, list) or isinstance(gene, np.ndarray):
        row = len(gene)
        fig = plt.figure(dpi=350,
                         constrained_layout=True,
                         figsize=(8, row * 2)
                         )

        gs = GridSpec(row, 1,
                      figure=fig)
        ax_list = []
        for index, value in enumerate(gene):
            ax = fig.add_subplot(gs[index, 0])
            freq_signal = adata[:, value].varm[domain]
            freq_signal = np.ravel(freq_signal)
            low = list(range(adata.uns['identify_svg_data']['fms_low'].shape[1]))
            plt.bar(low, freq_signal[low], color=colors[0])
            high = list(range(len(low), freq_signal.size))
            plt.bar(high, freq_signal[high], color=colors[1])
            ax.set_ylabel("signal")
            ax.spines['right'].set_color("none")
            ax.spines['top'].set_color("none")
            y_max = max(freq_signal)
            plt.ylim(0, y_max * 1.1)
            plt.xlim(0, freq_signal.size)
            plt.title("Gene: " + value)
            ax_list.append(ax)

        if save_path is not None:
            plt.savefig(f"{save_path}")
        plt.show()
        if return_fig:
            return ax_list


def ftu_freq_signal(adata,
                    ftu="ftu_1",
                    domain='freq_domain_svg',
                    figsize=(8, 2),
                    dpi=100,
                    color='#CA1C1C',
                    y_range=None,
                    xy_axis=True,
                    return_fig=False,
                    save_path=None,
                    **kwargs):
    if isinstance(ftu, str):
        # fig, ax = plt.subplots()
        plt.figure(figsize=figsize, dpi=dpi)
        ax = plt.gca()
        freq_signal = \
            adata.uns['detect_ftu_data']['freq_signal_ftu'].loc[ftu, :].values

        plt.bar(range(freq_signal.size), freq_signal, color=color)
        plt.grid(False)

        plt.title(ftu)
        ax.set_ylabel("signal")
        ax.spines['right'].set_color("none")
        ax.spines['top'].set_color("none")
        if y_range is not None:
            plt.ylim(y_range[0], y_range[1])
        else:
            plt.ylim(0, max(freq_signal) * 1.1)
        plt.xlim(0, freq_signal.size)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.set_ylabel("signal")
        if not xy_axis:
            ax.spines['left'].set_visible(False)
            ax.spines['bottom'].set_visible(False)
            ax.get_xaxis().set_visible(False)
            ax.get_yaxis().set_visible(False)
        if save_path is not None:
            plt.savefig(f"{save_path}")
        plt.show()
        if return_fig:
            return ax

    elif isinstance(ftu, list) or isinstance(ftu, np.ndarray):
        row = len(ftu)
        fig = plt.figure(dpi=350,
                         constrained_layout=True,
                         figsize=(8, row * 2)
                         )

        gs = GridSpec(row, 1,
                      figure=fig)
        ax_list = []
        for index, value in enumerate(ftu):
            ax = fig.add_subplot(gs[index, 0])
            freq_signal = \
                adata.uns['detect_ftu_data']['freq_signal_ftu'].loc[value,
                :].values

            plt.bar(range(freq_signal.size), freq_signal, color=color)
            plt.grid(False)

            plt.title(value)
            ax.set_ylabel("signal")
            ax.spines['right'].set_color("none")
            ax.spines['top'].set_color("none")
            if y_range != None:
                plt.ylim(y_range[0], y_range[1])
            else:
                plt.ylim(0, max(freq_signal) * 1.1)
            plt.xlim(0, freq_signal.size)
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)
            ax.set_ylabel("signal")
            if not xy_axis:
                ax.spines['left'].set_visible(False)
                ax.spines['bottom'].set_visible(False)
                ax.get_xaxis().set_visible(False)
                ax.get_yaxis().set_visible(False)
            ax_list.append(ax)

        if save_path is not None:
            plt.savefig(f"{save_path}")
        plt.show()
        if return_fig:
            return ax_list


def _subftu_freq_signal(adata,
                        subftu,
                        figsize=(6, 2),
                        dpi=100,
                        color='#CA1C1C',
                        y_range=[0, 0.1],
                        return_fig=False, **kwargs):
    # Show the frequency signal
    freq_signal = adata.uns['freq_signal_subftu'].loc[subftu, :].values
    plt.figure(figsize=figsize, dpi=dpi)
    low = list(range(len(freq_signal)))
    plt.bar(low, freq_signal, color=color)
    ax = plt.gca()
    plt.grid(False)
    ax.set_ylabel("signal")
    ax.spines['right'].set_color("none")
    ax.spines['top'].set_color("none")
    plt.ylim(y_range[0], y_range[1])
    plt.xlim(0, freq_signal.size)
    plt.title(subftu)
    plt.show()
    if return_fig:
        return ax


def gene_signal_umap(adata,
                     svg_list,
                     colors=['#C71C1C', '#BEBEBE'],
                     n_neighbors=15,
                     random_state=0,
                     return_fig=False,
                     save_path=None,
                     **kwargs):
    """
    The UMAP of svgs and non-svgs.

    Parameters
    ----------
    adata : anndata
        The spatial dataset.
    svg_list : list
        The svg list.
    colors : list, optional
        The colors corresponding to svgs and non-svgs. 
        The default is ['#C71C1C', '#BEBEBE'].
    n_neighbors : int, optional
        The neighbors when construct gene graph for umap.
        The default is 15.
    random_state : int, optional
        The ramdom state. The default is 0.
    return_fig : bool, optional
        Whether you need to return figure. The default is False.
    save_path : system path | None, optional
        The path including filename to save the figure.

    Returns
    -------
    fig : matploblib axes
        The figure.

    """
    low_length = adata.uns['identify_svg_data']['frequencies_low'].shape[0]
    freq_domain = adata.varm['freq_domain_svg'][:, :low_length].copy()
    freq_domain = preprocessing.normalize(freq_domain, norm='l1')
    freq_domain = pd.DataFrame(freq_domain)
    freq_domain.index = adata.var_names
    umap_adata = sc.AnnData(freq_domain)
    sc.pp.neighbors(umap_adata, n_neighbors=n_neighbors, use_rep='X')
    sc.tl.umap(umap_adata, random_state=0)
    adata.varm['freq_umap_svg'] = umap_adata.obsm['X_umap']
    print("""The umap coordinates of genes when identify svgs could be found in 
          adata.varm['freq_umap_svg']""")
    # svg_list
    umap_adata.obs['SpaGFT'] = 'Non-svgs'
    umap_adata.obs.loc[svg_list, 'SpaGFT'] = 'svgs'
    umap_adata.obs['SpaGFT'] = pd.Categorical(umap_adata.obs['SpaGFT'],
                                              categories=['svgs', 'Non-svgs'],
                                              ordered=True)
    umap_adata.uns['SpaGFT_colors'] = colors
    fig = sc.pl.umap(umap_adata, color="SpaGFT", return_fig=return_fig,
                     **kwargs)

    if save_path is not None:
        plt.savefig(f"{save_path}")
    plt.show()
    if return_fig:
        return fig


def _scatter_gene_distri(adata,
                         gene,
                         size=3,
                         shape='h',
                         cmap='magma',
                         spatial_info=['array_row', 'array_col'],
                         coord_ratio=0.7,
                         return_plot=False):
    if gene in adata.obs.columns:
        if isinstance(gene, str):
            plot_df = pd.DataFrame(adata.obs.loc[:, gene].values,
                                   index=adata.obs_names,
                                   columns=[gene])
        else:
            plot_df = pd.DataFrame(adata.obs.loc[:, gene],
                                   index=adata.obs_names,
                                   columns=gene)
        if spatial_info in adata.obsm_keys():
            plot_df['x'] = adata.obsm[spatial_info][:, 0]
            plot_df['y'] = adata.obsm[spatial_info][:, 1]
        elif set(spatial_info) <= set(adata.obs.columns):
            plot_coor = adata.obs
            plot_df['x'] = plot_coor.loc[:, spatial_info[0]].values
            plot_df['y'] = plot_coor.loc[:, spatial_info[1]].values

        if isinstance(gene, str):
            base_plot = (ggplot() + geom_point(plot_df, aes(x='x', y='y',
                                                            fill=gene),
                                               shape=shape, stroke=0.1,
                                               size=size) +
                         xlim(min(plot_df.x) - 1, max(plot_df.x) + 1) +
                         ylim(min(plot_df.y) - 1, max(plot_df.y) + 1) +
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
            for i in gene:
                base_plot = (ggplot() + geom_point(plot_df, aes(x='x', y='y',
                                                                fill=gene),
                                                   shape=shape, stroke=0.1,
                                                   size=size) +
                             xlim(min(plot_df.x) - 1, max(plot_df.x) + 1) +
                             ylim(min(plot_df.y) - 1, max(plot_df.y) + 1) +
                             scale_fill_cmap(cmap_name=cmap) +
                             coord_equal(ratio=coord_ratio) +
                             theme_classic() +
                             theme(legend_position=('right'),
                                   legend_background=element_blank(),
                                   legend_key_width=4,
                                   legend_key_height=50)
                             )
                print(base_plot)

        return
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
    if isinstance(gene, str):
        base_plot = (ggplot() + geom_point(plot_df, aes(x='x', y='y',
                                                        fill=gene),
                                           shape=shape, stroke=0.1,
                                           size=size) +
                     xlim(min(plot_df.x) - 1, max(plot_df.x) + 1) +
                     ylim(min(plot_df.y) - 1, max(plot_df.y) + 1) +
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
        for i in gene:
            base_plot = (ggplot() + geom_point(plot_df, aes(x='x', y='y',
                                                            fill=gene),
                                               shape=shape, stroke=0.1,
                                               size=size) +
                         xlim(min(plot_df.x) - 1, max(plot_df.x) + 1) +
                         ylim(min(plot_df.y) - 1, max(plot_df.y) + 1) +
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


def _umap_svg_cluster(adata,
                      svg_list=None,
                      size=None,
                      coord_ratio=1,
                      return_plot=True):
    if svg_list == None:
        ftup_df = adata.var.copy()
        svg_list = ftup_df[ftup_df.cutoff_gft_score][ftup_df.qvalue < 0.05].index
    plot_df = adata.uns['gft_umap_ftu']
    plot_df = pd.DataFrame(plot_df)
    plot_df.index = adata.var_names
    plot_df.columns = ['UMAP_1', 'UMAP_2']

    plot_df.loc[svg_list, 'gene'] = 'svg'
    plot_df['gene'] = pd.Categorical(plot_df['gene'],
                                     categories=['svg', 'Non-svg'],
                                     ordered=True)
    plot_df['radius'] = size
    # plot
    if size is None:
        base_plot = (ggplot(plot_df, aes(x='UMAP_1', y='UMAP_2', fill='gene')) +
                     geom_point(color='white', stroke=0.25) +
                     scale_fill_manual(values=colors) +
                     theme_classic() +
                     coord_equal(ratio=coord_ratio))
    else:
        base_plot = (ggplot(plot_df, aes(x='UMAP_1', y='UMAP_2', fill='gene')) +
                     geom_point(size=size, color='white', stroke=0.25) +
                     scale_fill_manual(values=colors) +
                     theme_classic() +
                     coord_equal(ratio=coord_ratio))
    print(base_plot)
    if return_plot:
        return base_plot


def _scatter_ftu_binary(adata,
                        ftu,
                        size=3,
                        shape='h',
                        spatial_info=['array_row', 'array_col'],
                        colors=['#CA1C1C', '#CCCCCC'],
                        coord_ratio=0.7,
                        return_plot=False):
    if '-' in ftu:
        ftu = 'ftu-' + ftu.split('-')[0] + "_subftu-" + ftu.split('-')[1]
        plot_df = adata.obsm['subftu_binary']
    else:
        ftu = 'ftu_' + ftu
        plot_df = adata.obsm['ftu_binary']
    plot_df = pd.DataFrame(plot_df)
    if spatial_info in adata.obsm_keys():
        plot_df['x'] = adata.obsm[spatial_info][:, 0]
        plot_df['y'] = adata.obsm[spatial_info][:, 1]
    elif set(spatial_info) <= set(adata.obs.columns):
        plot_coor = adata.obs
        plot_df = plot_df[ftu]
        plot_df = pd.DataFrame(plot_df)
        plot_df['x'] = plot_coor.loc[:, spatial_info[0]].values
        plot_df['y'] = plot_coor.loc[:, spatial_info[1]].values
    plot_df['radius'] = size
    plot_df[ftu] = plot_df[ftu].values.astype(int)
    plot_df[ftu] = plot_df[ftu].values.astype(str)
    plot_df[ftu] = pd.Categorical(plot_df[ftu],
                                  categories=['1', '0'],
                                  ordered=True)
    base_plot = (ggplot() + geom_point(plot_df, aes(x='x', y='y', fill=ftu),
                                       shape=shape, stroke=0.1, size=size) +
                 xlim(min(plot_df.x) - 1, max(plot_df.x) + 1) +
                 ylim(min(plot_df.y) - 1, max(plot_df.y) + 1) +
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


def umap_svg(adata,
             svg_list=None,
             colors=['#CA1C1C', '#CCCCCC'],
             size=2,
             coord_ratio=0.7,
             return_plot=False):
    if 'gft_umap_svg' not in adata.varm_keys():
        raise KeyError("Please run SpaGFT.calculate_frequency_domain firstly")
    plot_df = adata.varm['gft_umap_svg']
    plot_df = pd.DataFrame(plot_df)
    plot_df.index = adata.var_names
    plot_df.columns = ['UMAP_1', 'UMAP_2']
    plot_df['gene'] = 'Non-svg'
    if svg_list == None:
        ftup_df = adata.var.copy()
        svg_list = ftup_df[ftup_df.cutoff_gft_score][ftup_df.qvalue < 0.05].index
    plot_df.loc[svg_list, 'gene'] = 'svg'
    plot_df['gene'] = pd.Categorical(plot_df['gene'],
                                     categories=['svg', 'Non-svg'],
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


def _visualize_fms(adata,
                   rank=1,
                   low=True,
                   size=3,
                   cmap='magma',
                   spatial_info=['array_row', 'array_col'],
                   shape='h',
                   coord_ratio=1,
                   return_plot=False):
    if low:
        plot_df = pd.DataFrame(adata.uns['fms_low'])
        plot_df.index = adata.obs.index
        plot_df.columns = ['low_FM_' + \
                           str(i + 1) for i in range(plot_df.shape[1])]
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
                                           shape=shape,
                                           stroke=0.1,
                                           size=size) +
                     xlim(min(plot_df.x) - 1, max(plot_df.x) + 1) +
                     ylim(min(plot_df.y) - 1, max(plot_df.y) + 1) +
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
        plot_df.columns = ['high_FM_' + str(i + 1) for i in \
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
                                                             str(adata.uns['fms_high'].shape[1] - \
                                                                 rank + 1)),
                                           shape=shape, stroke=0.1, size=size) +
                     xlim(min(plot_df.x) - 1, max(plot_df.x) + 1) +
                     ylim(min(plot_df.y) - 1, max(plot_df.y) + 1) +
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


def ftu_heatmap_signal_ftu_id_card(adata,
                                   ftu,
                                   ax=None,
                                   title=None):
    freq_signal = \
        adata.uns['detect_ftu_data']['freq_signal_ftu'].loc[ftu,
        :].values.reshape(1, -1)
    # print(freq_signal)
    if title != None:
        plt.title(title, fontsize=10)
    sns.heatmap(data=freq_signal, cbar=False, cmap="Reds")
    ax.tick_params(left=False, bottom=False)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.get_xaxis().set_visible(False)
    ax.get_yaxis().set_visible(False)


# ftu id card
def cell_type_proportion_box_ftu_id_card(cell_type_name,
                                         cell_type_proportion_data,
                                         ax=None,
                                         title=None):
    boxplot_data = []
    for val in cell_type_name:
        boxplot_data.append(cell_type_proportion_data[val])
    labels = [x.replace("q05cell_abundance_w_sf_", "") for x in cell_type_name]
    plt.title(title, fontsize=20)
    ax.boxplot(boxplot_data, labels=labels, showfliers=False,
               patch_artist=True,
               boxprops={"facecolor": "#FF2A6B"},
               medianprops={"color": "black"})
    ax.yaxis.set_tick_params(labelsize=7)
    ax.xaxis.set_tick_params(labelsize=7)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    for tick in ax.get_xticklabels():
        tick.set_rotation(90)


# ftu id card
def ftu_spatial_map_scatter_ftu_id_card(adata,
                                        ftu_name,
                                        ftu_color,
                                        title,
                                        shape='h',
                                        radius=0.5,
                                        spatial_info=['array_row', 'array_col'],
                                        ax=None):
    x = []
    y = []
    ftu = list(adata.obsm["ftu_binary"][ftu_name].values)
    ftu = [int(x) for x in ftu]
    cmap = ListedColormap(["lightgray", ftu_color])

    if spatial_info in adata.obsm_keys():
        x = adata.obsm[spatial_info][:, 1]
        y = adata.obsm[spatial_info][:, 0]
    elif set(spatial_info) <= set(adata.obs.columns):
        plot_coor = adata.obs
        x = plot_coor.loc[:, spatial_info[0]].values
        y = plot_coor.loc[:, spatial_info[1]].values

    if title != None:
        plt.title(title, y=-0.2, fontsize=10)
    ax.scatter(y, max(x) - x, s=radius, c=ftu, cmap=cmap, marker=shape)
    ax.minorticks_on()
    ax.yaxis.set_tick_params(labelsize=10)
    ax.xaxis.set_tick_params(labelsize=10)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.get_xaxis().set_visible(False)
    ax.get_yaxis().set_visible(False)
    ax.set_aspect('equal')
    ax.grid(False)


# ftu id card
def ftu_overlapped_scatter_ftu_id_card(adata,
                                       ftu_1_code,
                                       ftu_2_code,
                                       ftu_1_color,
                                       ftu_2_color,
                                       overlapped_color,
                                       title,
                                       marker='h',
                                       radius=0.03,
                                       spatial_info=['array_row', 'array_col'],
                                       ax=None):
    ftu_1 = adata.obsm["ftu_binary"][f"ftu_{ftu_1_code}"].values
    ftu_2 = adata.obsm["ftu_binary"][f"ftu_{ftu_2_code}"].values
    x = []
    y = []

    if spatial_info in adata.obsm_keys():
        x = adata.obsm[spatial_info][:, 1]
        y = adata.obsm[spatial_info][:, 0]
    elif set(spatial_info) <= set(adata.obs.columns):
        plot_coor = adata.obs
        x = plot_coor.loc[:, spatial_info[0]].values
        y = plot_coor.loc[:, spatial_info[1]].values
    plt.title(title, y=-0.2)

    ftu_x = max(x) - x
    ftu_1_value = []
    ftu_2_value = []
    ftu_overlapped_value = []
    ftu_blank_value = []
    for index in range(len(ftu_1)):
        if ftu_1[index] == ftu_2[index] and ftu_1[index] == "1":
            ftu_overlapped_value.append([y[index], ftu_x[index], 3])
        elif ftu_1[index] == "1":
            ftu_1_value.append([y[index], ftu_x[index], 1])
        elif ftu_2[index] == "1":
            ftu_2_value.append([y[index], ftu_x[index], 2])
        else:
            ftu_blank_value.append([y[index], ftu_x[index], 0])
    label_name = [f"ftu_{ftu_1_code}", f"ftu_{ftu_2_code}", "ftu_overlap", "ftu_blank"]
    ftu_1_value = np.array(ftu_1_value)
    ftu_overlapped_value = np.array(ftu_overlapped_value)
    ftu_2_value = np.array(ftu_2_value)
    ftu_blank_value = np.array(ftu_blank_value)
    scatter_1 = ax.scatter(ftu_1_value[:, 0], ftu_1_value[:, 1],
                           marker=marker, s=radius, c=ftu_1_color)
    scatter_2 = ax.scatter(ftu_2_value[:, 0], ftu_2_value[:, 1],
                           marker=marker, s=radius, c=ftu_2_color)
    scatter_overlap = ax.scatter(ftu_overlapped_value[:, 0],
                                 ftu_overlapped_value[:, 1],
                                 marker=marker,
                                 s=radius,
                                 c=overlapped_color)
    scatter_blank = ax.scatter(ftu_blank_value[:, 0],
                               ftu_blank_value[:, 1],
                               marker=marker,
                               s=radius,
                               c="lightgray")
    ax.minorticks_on()
    ax.yaxis.set_tick_params(labelsize=10)
    ax.xaxis.set_tick_params(labelsize=10)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.get_xaxis().set_visible(False)
    ax.get_yaxis().set_visible(False)
    ax.set_aspect('equal')
    ax.legend([scatter_1], [f"ftu_{ftu_1_code}"])
    ax.legend([scatter_2], [f"ftu_{ftu_2_code}"])
    ax.legend([scatter_overlap], ["ftu_overlap"])
    ax.legend([scatter_blank], ["ftu_blank"])

    plt.legend(label_name, fontsize=5, loc=7, frameon=False,
               bbox_to_anchor=(1, 0, 0.35, 1))
    ax.grid(False)


# ftu id card
def scatter_svgs_distri_ftu_id_card(adata,
                                    gene,
                                    size=None,
                                    shape='h',
                                    cmap='magma',
                                    spatial_info=['array_row', 'array_col'],
                                    ax=None,
                                    coord_ratio=1,
                                    return_plot=False):
    if gene in adata.obs.columns:
        plot_df = pd.DataFrame(adata.obs.loc[:, gene].values,
                               index=adata.obs_names,
                               columns=[gene])
        if spatial_info in adata.obsm_keys():
            plot_df['x'] = adata.obsm[spatial_info][:, 1]
            plot_df['y'] = adata.obsm[spatial_info][:, 0]
        elif set(spatial_info) <= set(adata.obs.columns):
            plot_coor = adata.obs
            plot_df = plot_df[gene]
            plot_df = pd.DataFrame(plot_df)
            plot_df['x'] = plot_coor.loc[:, spatial_info[0]].values
            plot_df['y'] = plot_coor.loc[:, spatial_info[1]].values
        plot_df['radius'] = size
        if size is not None:
            plot_scatter_ftu_id_card(plot_df.y,
                                     max(plot_df.x) - plot_df.x, plot_df[gene],
                                     gene, cmap, radius=plot_df['radius'], ax=ax)
        else:
            plot_scatter_ftu_id_card(plot_df.y,
                                     max(plot_df.x) - plot_df.x, plot_df[gene],
                                     gene, cmap, ax=ax)
        return
    if ss.issparse(adata.X):
        plot_df = pd.DataFrame(adata.X.todense(), index=adata.obs_names,
                               columns=adata.var_names)
    else:
        plot_df = pd.DataFrame(adata.X, index=adata.obs_names,
                               columns=adata.var_names)
    if spatial_info in adata.obsm_keys():
        plot_df['x'] = adata.obsm[spatial_info][:, 1]
        plot_df['y'] = adata.obsm[spatial_info][:, 0]
    elif set(spatial_info) <= set(adata.obs.columns):
        plot_coor = adata.obs
        plot_df = plot_df[gene]
        plot_df = pd.DataFrame(plot_df)
        plot_df['x'] = plot_coor.loc[:, spatial_info[0]].values
        plot_df['y'] = plot_coor.loc[:, spatial_info[1]].values
    plot_df['radius'] = size
    plot_df = plot_df.sort_values(by=gene, ascending=True)
    # print(plot_df)
    if size is not None:
        plot_scatter_ftu_id_card(plot_df.y, max(plot_df.x) - plot_df.x,
                                 plot_df[gene], gene,
                                 cmap, radius=plot_df['radius'], ax=ax)
    else:
        plot_scatter_ftu_id_card(plot_df.y, max(plot_df.x) - plot_df.x,
                                 plot_df[gene], gene,
                                 cmap, ax=ax)


# ftu id card
def ftu_freq_signal_ftu_id_card(adata,
                                ftu,
                                color='#CA1C1C',
                                y_range=[0, 0.08],
                                ax=None, title=None):
    # Show the frequency signal
    freq_signal = \
        adata.uns['detect_ftu_data']['freq_signal_ftu'].loc[ftu, :].values
    low = list(range(len(freq_signal)))
    plt.bar(low, freq_signal, color=color)
    plt.grid(False)
    if title != None:
        plt.title(title, fontsize=20)
    ax.spines['right'].set_color("none")
    ax.spines['top'].set_color("none")
    plt.ylim(y_range[0], y_range[1])
    plt.xlim(0, freq_signal.size)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.get_xaxis().set_visible(False)
    ax.get_yaxis().set_visible(False)


def plot_scatter_ftu_id_card(x, y, colors, title, cmap, shape='h',
                             radius=None, ax=None, up_title=False):
    # fig, ax = plt.subplots()
    # fig.subplots_adjust(right=0.9)
    if up_title:
        plt.title(title)
    else:
        plt.title(title, y=-0.2)
    if radius is not None:
        scatter = ax.scatter(x, y, s=radius, c=colors, cmap=cmap, marker=shape)
    else:
        scatter = ax.scatter(x, y, c=colors, cmap=cmap, marker=shape)
    ax.minorticks_on()
    ax.yaxis.set_tick_params(labelsize=10)
    ax.xaxis.set_tick_params(labelsize=10)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.get_xaxis().set_visible(False)
    ax.get_yaxis().set_visible(False)
    ax.set_aspect('equal')
    ax.grid(False)
    return scatter


def scatter_gene(adata,
                 gene,
                 size=None,
                 shape='o',
                 cmap='magma',
                 spatial_info=['array_row', 'array_col'],
                 dpi=100,
                 return_fig=False,
                 save_path=None):
    """
    Visualize genes through given spatial information.

    Parameters
    ----------
    adata : anndata
        spatial datasets. 
    gene : str or list
        The genes which will be visualized.
    size : float, optional
        The size of spots. 
        The default is None.
    shape : str, optional
        The shape of the spots. 
        The default is 'o'.
    cmap : str, optional
        The color theme used. 
        The default is "magma".
    spatial_info : str or list, optional
        The spatial information key in adata.obsm or columns in adata.obs. 
        The default is ['array_row', 'array_col'].
    dpi : int, optional
        DPI. The default is 100.
    return_fig : bool, optional
        Where to return the figure.
        The default is False.
    save_path : system path, optional
        The path for the saving figure.
        The default is None.
    Raises
    ------
    KeyError
        gene should be adata.var_names.

    """
    if isinstance(gene, np.ndarray):
        gene = list(gene)
    if isinstance(gene, pd.core.indexes.base.Index):
        gene = list(gene)
    if ss.issparse(adata.X):
        if isinstance(gene, str):
            plot_df = pd.DataFrame(adata[:, gene].X.todense(),
                                   index=adata.obs_names,
                                   columns=[gene])
        elif isinstance(gene, list) or isinstance(gene, np.ndarray):
            plot_df = pd.DataFrame(adata[:, gene].X.todense(),
                                   index=adata.obs_names,
                                   columns=gene)
        else:
            raise KeyError(f"{gene} is invalid!")
    else:
        if isinstance(gene, str):
            plot_df = pd.DataFrame(adata[:, gene].X,
                                   index=adata.obs_names,
                                   columns=[gene])
        elif isinstance(gene, list) or isinstance(gene, np.ndarray):
            plot_df = pd.DataFrame(adata[:, gene].X,
                                   index=adata.obs_names,
                                   columns=gene)
        else:
            raise KeyError(f"{gene} is invalid!")
    if spatial_info in adata.obsm_keys():
        plot_df['x'] = adata.obsm[spatial_info][:, 1]
        plot_df['y'] = adata.obsm[spatial_info][:, 0]
    elif set(spatial_info) <= set(adata.obs.columns):
        plot_coor = adata.obs
        plot_df = plot_df[gene]
        plot_df = pd.DataFrame(plot_df)
        plot_df['x'] = plot_coor.loc[:, spatial_info[0]].values
        plot_df['y'] = plot_coor.loc[:, spatial_info[1]].values
        # print(plot_df)
    if isinstance(gene, str):
        fig, ax = plt.subplots()
        if size == None:
            scatter = plot_scatter_ftu_id_card(x=plot_df.y,
                                               y=max(plot_df.x) - plot_df.x,
                                               colors=plot_df[gene],
                                               title=gene,
                                               shape=shape,
                                               cmap=cmap,
                                               ax=ax,
                                               up_title=True)
            plt.colorbar(scatter, ax=ax)
        elif isinstance(size, int) or isinstance(size, float):
            scatter = plot_scatter_ftu_id_card(x=plot_df.y,
                                               y=max(plot_df.x) - plot_df.x,
                                               colors=plot_df[gene],
                                               title=gene,
                                               shape=shape,
                                               cmap=cmap,
                                               radius=size,
                                               ax=ax,
                                               up_title=True)
            plt.colorbar(scatter, ax=ax)

        if save_path is not None:
            plt.savefig(save_path)
        plt.show()
        if return_fig:
            return ax

    elif isinstance(gene, list) or isinstance(gene, np.ndarray):
        row = math.ceil(len(gene) / 4)
        fig = plt.figure(dpi=dpi,
                         constrained_layout=True,
                         figsize=(20, row * 5))

        gs = GridSpec(row, 4,
                      figure=fig)
        ax_list = []
        for index, value in enumerate(gene):
            ax = fig.add_subplot(gs[index // 4, index % 4])

            if size == None:
                scatter = plot_scatter_ftu_id_card(x=plot_df.y,
                                                   y=max(plot_df.x) - plot_df.x,
                                                   colors=plot_df[value],
                                                   title=value,
                                                   shape=shape,
                                                   cmap=cmap,
                                                   ax=ax,
                                                   up_title=True)
                plt.colorbar(scatter, ax=ax)
            elif isinstance(size, int) or isinstance(size, float):
                scatter = plot_scatter_ftu_id_card(x=plot_df.y,
                                                   y=max(plot_df.x) - plot_df.x,
                                                   colors=plot_df[value],
                                                   title=value,
                                                   shape=shape,
                                                   cmap=cmap,
                                                   radius=size,
                                                   ax=ax,
                                                   up_title=True)
                plt.colorbar(scatter, ax=ax)
            ax_list.append(ax)

        if save_path:
            plt.savefig(save_path)
        plt.show()
        if return_fig:
            return ax_list


def scatter_ftu(adata,
                ftu,
                shape='o',
                ftu_color='#FF6879',
                size=None,
                spatial_info=['array_row', 'array_col'],
                save_path=None,
                return_fig=False):
    """
    Plot the spatial map of a ftu.

    Parameters
    ----------
    adata : anndata
        spatial datasets. 
    ftu : str
        The ftu indicator.
    shape : str, optional
        The shape of spots. The default is 'o'.
    ftu_color : str, optional
        The color of the ftu(s). The default is '#FF6879'.
    size : float, optional
        The size of spots. The default is None.
    spatial_info : str or list, optional
        The spatial information key in adata.obsm or columns in adata.obs. 
        The default is ['array_row', 'array_col'].
    return_fig : bool, optional
        Where to return the figure. The default is False.
    save_path : system path, optional
        The path for the saving figure. The default is None.

    Returns
    -------
    ax_list 
        If return_fig == True, the figure objects will be returned.

    """
    x = []
    y = []
    if spatial_info in adata.obsm_keys():
        x = adata.obsm[spatial_info][:, 1]
        y = adata.obsm[spatial_info][:, 0]
    elif set(spatial_info) <= set(adata.obs.columns):
        plot_coor = adata.obs
        x = plot_coor.loc[:, spatial_info[0]].values
        y = plot_coor.loc[:, spatial_info[1]].values
    if isinstance(ftu, str):
        ftu_value = [int(x) for x in list(adata.obsm["ftu_binary"][ftu].values)]
        fig, ax = plt.subplots()
        cmap_ftu = ListedColormap(["#b4b4b4", ftu_color])

        plt.title(ftu)
        if size is not None:
            scatter = ax.scatter(y, max(x) - x, s=size, c=ftu_value,
                                 cmap=cmap_ftu, marker=shape)
        else:
            scatter = ax.scatter(y, max(x) - x, c=ftu_value,
                                 cmap=cmap_ftu, marker=shape)
        ax.minorticks_on()
        ax.yaxis.set_tick_params(labelsize=10)
        ax.xaxis.set_tick_params(labelsize=10)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['left'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
        ax.get_xaxis().set_visible(False)
        ax.get_yaxis().set_visible(False)
        ax.set_aspect('equal')
        ax.grid(False)
        plt.legend(*scatter.legend_elements(), loc="center right",
                   bbox_to_anchor=(1, 0, 0.15, 1))

        if save_path:
            plt.savefig(f"{save_path}")
        plt.show()
        if return_fig:
            return ax
    elif isinstance(ftu, list) or isinstance(ftu, np.ndarray):
        row = math.ceil(len(ftu) / 4)
        fig = plt.figure(dpi=200,
                         constrained_layout=True,
                         figsize=(20, row * 5)
                         )

        gs = GridSpec(row,
                      4,
                      figure=fig)
        # print(row, 4)
        ax_list = []
        for index, value in enumerate(ftu):
            ax = fig.add_subplot(gs[index // 4, index % 4])
            ftu_value = [int(x) for x in \
                         list(adata.obsm["ftu_binary"][value].values)]
            cmap_ftu = ListedColormap(["#b4b4b4", ftu_color])

            plt.title(value)
            if size is not None:
                scatter = ax.scatter(y, max(x) - x, s=size, c=ftu_value,
                                     cmap=cmap_ftu, marker=shape)
            else:
                scatter = ax.scatter(y, max(x) - x, c=ftu_value,
                                     cmap=cmap_ftu, marker=shape)

            ax.minorticks_on()
            ax.yaxis.set_tick_params(labelsize=10)
            ax.xaxis.set_tick_params(labelsize=10)
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)
            ax.spines['left'].set_visible(False)
            ax.spines['bottom'].set_visible(False)
            ax.get_xaxis().set_visible(False)
            ax.get_yaxis().set_visible(False)
            ax.set_aspect('equal')
            ax.grid(False)
            plt.legend(*scatter.legend_elements(), loc="center right",
                       bbox_to_anchor=(1, 0, 0.15, 1))
            ax_list.append(ax)

        if save_path:
            plt.savefig(f"{save_path}")
        plt.show()
        if return_fig:
            return ax_list


def scatter_ftu_expression(adata,
                           ftu,
                           cmap='magma',
                           shape='o',
                           size=None,
                           spatial_info=['array_row', 'array_col'],
                           save_path=None,
                           return_fig=False):
    x = []
    y = []
    if spatial_info in adata.obsm_keys():
        x = adata.obsm[spatial_info][:, 1]
        y = adata.obsm[spatial_info][:, 0]
    elif set(spatial_info) <= set(adata.obs.columns):
        plot_coor = adata.obs
        x = plot_coor.loc[:, spatial_info[0]].values
        y = plot_coor.loc[:, spatial_info[1]].values
    if isinstance(ftu, str):
        ftu_value = adata.obsm["ftu_pseudo_expression"][ftu].values
        fig, ax = plt.subplots()
        plt.title(ftu)
        if size is not None:
            scatter = ax.scatter(y, max(x) - x, s=size, c=ftu_value, cmap=cmap,
                                 marker=shape)
        else:
            scatter = ax.scatter(y, max(x) - x, c=ftu_value, cmap=cmap,
                                 marker=shape)
        ax.minorticks_on()
        ax.yaxis.set_tick_params(labelsize=10)
        ax.xaxis.set_tick_params(labelsize=10)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['left'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
        ax.get_xaxis().set_visible(False)
        ax.get_yaxis().set_visible(False)
        ax.set_aspect('equal')
        plt.colorbar(scatter, ax=ax)
        ax.grid(False)

        if save_path:
            plt.savefig(f"{save_path}")
        plt.show()
        if return_fig:
            return ax
    elif isinstance(ftu, list) or isinstance(ftu, np.ndarray):
        row = math.ceil(len(ftu) / 4)
        fig = plt.figure(dpi=350,
                         constrained_layout=True,
                         figsize=(20, row * 5)
                         )

        gs = GridSpec(row, 4,
                      figure=fig)
        # print(row, 4)
        ax_list = []
        for index, value in enumerate(ftu):
            ax = fig.add_subplot(gs[index // 4, index % 4])
            ftu_value = adata.obsm["ftu_pseudo_expression"][value].values

            plt.title(value)
            if size is not None:
                scatter = ax.scatter(y, max(x) - x, s=size, c=ftu_value,
                                     cmap=cmap, marker=shape)
            else:
                scatter = ax.scatter(y, max(x) - x, c=ftu_value, cmap=cmap,
                                     marker=shape)
            ax.minorticks_on()
            ax.yaxis.set_tick_params(labelsize=10)
            ax.xaxis.set_tick_params(labelsize=10)
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)
            ax.spines['left'].set_visible(False)
            ax.spines['bottom'].set_visible(False)
            ax.get_xaxis().set_visible(False)
            ax.get_yaxis().set_visible(False)
            ax.set_aspect('equal')
            plt.colorbar(scatter, ax=ax)
            ax.grid(False)
            ax_list.append(ax)

        if save_path:
            plt.savefig(f"{save_path}")
        plt.show()
        if return_fig:
            return ax_list


def overlap_curve(adata, save_path=None, return_fig=False):
    curve_dot = adata.uns['detect_ftu_data']['overlap_curve']
    curve_x = [float(x.replace("res_", "")) for x in list(curve_dot.index)]
    curve_y = list(curve_dot["score"].values)
    x_values = curve_x  # .tolist()
    y_values = curve_y
    ax = plt.scatter(x_values, y_values)
    plt.plot(x_values, y_values)
    plt.title("Overlap curve", fontsize=24)
    plt.xlabel("Resolution", fontsize=14)
    plt.ylabel("Score", fontsize=14)

    if save_path is not None:
        plt.savefig(f"{save_path}")
    plt.show()
    if return_fig:
        return ax


def scatter_umap_clustering(adata,
                            svg_list,
                            size=3,
                            alpha=1,
                            return_fig=False,
                            save_path=None):
    """
    Visualize clustering results.

    Parameters
    ----------
    adata : anndata
        spatial datasets. 
    svg_list : list
        svg list.
    size : float
        The size of genes.
        The default is 3.
    alpha : float, optional
        Transparency.
        The default is 1.
    return_fig : bool, optional
        Where to return the figure. The default is False.
    save_path : system path, optional
        The path for the saving figure. The default is None.

    Raises
    ------
    KeyError
        Run sgagft.identify_tm() before this step.

    Returns
    -------
    base_plot : plotnine object
        The figure.

    """
    current_genes = adata.uns['detect_ftu_data']['gft_umap_ftu'].index.tolist()
    if set(svg_list) <= set(current_genes):
        svg_list = np.intersect1d(svg_list, current_genes)
    else:
        diff_genes = np.setdiff1d(svg_list, current_genes)
        raise KeyError(f"{diff_genes} are not calculated in the above step.")
    plot_df = pd.concat((adata.uns['detect_ftu_data'] \
                             ['gft_umap_ftu'].loc[svg_list, :],
                         adata.var.loc[svg_list, :].ftu), axis=1)

    categories = [eval(i) for i in np.unique(plot_df.ftu)]
    categories = np.sort(np.array(categories))
    categories = categories.astype(str)
    plot_df.ftu = pd.Categorical(plot_df.ftu,
                                 categories=categories)
    base_plot = (ggplot(plot_df,
                        aes('UMAP_1', 'UMAP_2',
                            fill='ftu'))
                 + geom_point(size=size, alpha=alpha, stroke=0.1)
                 + scale_fill_hue(s=0.9, l=0.65, h=0.0417, color_space='husl')
                 + theme_classic())
    print(base_plot)

    if save_path is not None:
        base_plot.save(f"{save_path}")
    if return_fig:
        return base_plot


def scatter_ftu_gene(adata,
                     ftu,
                     gene,  # list
                     shape='o',
                     cmap="magma",
                     ftu_color='#FF6879',
                     size=None,
                     spatial_info=['array_row', 'array_col'],
                     return_fig=False,
                     save_path=None):
    """
    Plot a ftu and several genes simultaneously.

    Parameters
    ----------
    adata : anndata
        spatial datasets. 
    ftu : str
        The ftu indicator.
    gene : list
        The list of gene names.
    shape : str, optional
        The shape of spots.
        The default is 'o'.
    cmap : str, optional
        The color theme used. 
        The default is "magma".
    ftu_color : str, optional
        The color used for representing the ftu. 
        The default is '#FF6879'.
    size : float, optional
        The size of spots. The default is None.
    spatial_info : str or list, optional
        The spatial information key in adata.obsm or columns in adata.obs. 
        The default is ['array_row', 'array_col'].
    return_fig : bool, optional
        Where to return the figure. The default is False.
    save_path : system path, optional
        The path for the saving figure. The default is None.
    Raises
    ------
    KeyError
        gene should be found at adata.var_names.

    Returns
    -------
    ax_list : matplotlist subaxes object list
        The figures.

    """

    if isinstance(gene, str):
        gene = [gene]
    if isinstance(gene, pd.core.indexes.base.Index):
        gene = list(gene)
    x = []
    y = []
    if spatial_info in adata.obsm_keys():
        x = adata.obsm[spatial_info][:, 1]
        y = adata.obsm[spatial_info][:, 0]
    elif set(spatial_info) <= set(adata.obs.columns):
        plot_coor = adata.obs
        x = plot_coor.loc[:, spatial_info[0]].values
        y = plot_coor.loc[:, spatial_info[1]].values
    row = math.ceil((1 + len(gene)) / 4)
    fig = plt.figure(dpi=350,
                     constrained_layout=True,
                     figsize=(20, row * 5)
                     )

    gs = GridSpec(row, 4,
                  figure=fig)
    ax_list = []
    ###########################################################
    ax_ftu = fig.add_subplot(gs[0, 0])
    ftu_value = [int(x) for x in list(adata.obsm["ftu_binary"][ftu].values)]
    cmap_ftu = ListedColormap(["#b4b4b4", ftu_color])

    plt.title(ftu)
    if size is not None:
        scatter = ax_ftu.scatter(y, max(x) - x, s=size, c=ftu_value, cmap=cmap_ftu)
    else:
        scatter = ax_ftu.scatter(y, max(x) - x, c=ftu_value, cmap=cmap_ftu)

    ax_ftu.minorticks_on()
    ax_ftu.yaxis.set_tick_params(labelsize=10)
    ax_ftu.xaxis.set_tick_params(labelsize=10)
    ax_ftu.spines['top'].set_visible(False)
    ax_ftu.spines['right'].set_visible(False)
    ax_ftu.spines['left'].set_visible(False)
    ax_ftu.spines['bottom'].set_visible(False)
    ax_ftu.get_xaxis().set_visible(False)
    ax_ftu.get_yaxis().set_visible(False)
    ax_ftu.set_aspect('equal')
    ax_ftu.grid(False)
    plt.legend(*scatter.legend_elements(), loc="center right",
               bbox_to_anchor=(1, 0, 0.15, 1))
    ax_list.append(ax_ftu)
    #########################

    if isinstance(gene, np.ndarray):
        gene = list(gene)
    if ss.issparse(adata.X):
        if isinstance(gene, str):
            plot_df = pd.DataFrame(adata[:, gene].X.todense(),
                                   index=adata.obs_names,
                                   columns=[gene])
        elif isinstance(gene, list) or isinstance(gene, np.ndarray):
            plot_df = pd.DataFrame(adata[:, gene].X.todense(),
                                   index=adata.obs_names,
                                   columns=gene)
        else:
            raise KeyError(f"{gene} is invalid!")
    else:
        if isinstance(gene, str):
            plot_df = pd.DataFrame(adata[:, gene].X,
                                   index=adata.obs_names,
                                   columns=[gene])
        elif isinstance(gene, list) or isinstance(gene, np.ndarray):
            plot_df = pd.DataFrame(adata[:, gene].X,
                                   index=adata.obs_names,
                                   columns=gene)
        else:
            raise KeyError(f"{gene} is invalid!")
    if spatial_info in adata.obsm_keys():
        plot_df['x'] = adata.obsm[spatial_info][:, 1]
        plot_df['y'] = adata.obsm[spatial_info][:, 0]
    elif set(spatial_info) <= set(adata.obs.columns):
        plot_coor = adata.obs
        plot_df = plot_df[gene]
        plot_df = pd.DataFrame(plot_df)
        plot_df['x'] = plot_coor.loc[:, spatial_info[0]].values
        plot_df['y'] = plot_coor.loc[:, spatial_info[1]].values
        # print(plot_df)
    if isinstance(gene, list) or isinstance(gene, np.ndarray):
        for index, value in enumerate(gene):
            ax = fig.add_subplot(gs[(1 + index) // 4, (1 + index) % 4])

            if size == None:
                scatter = plot_scatter_ftu_id_card(x=plot_df.y,
                                                   y=max(plot_df.x) - plot_df.x,
                                                   shape=shape,
                                                   colors=plot_df[value],
                                                   title=value,
                                                   cmap=cmap,
                                                   ax=ax,
                                                   up_title=True)
            elif isinstance(size, int) or isinstance(size, float):
                scatter = plot_scatter_ftu_id_card(x=plot_df.y,
                                                   y=max(plot_df.x) - plot_df.x,
                                                   colors=plot_df[value],
                                                   title=value,
                                                   cmap=cmap,
                                                   radius=size,
                                                   ax=ax,
                                                   up_title=True)
            plt.colorbar(scatter, ax=ax)
            ax_list.append(ax)

    if save_path:
        plt.savefig(save_path)
    plt.show()
    if return_fig:
        return ax_list


def draw_ftu_id_card(adata,
                     svg_list,
                     ftu,
                     deconvolution_key='cell_type_proportion',
                     spatial_info=['array_row', 'array_col'],
                     shape='h',
                     dpi=350,
                     size=[7, 0.8],
                     return_fig=False,
                     save_path=None):
    """
    Plot the details of a ftu to generate a ftu ID CARD.

    Parameters
    ----------
    adata : anndata
        spatial datasets. 
    svg_list : list
        svg list.
    ftu : str
        The ftu indicator.
    deconvolution_key : str, optional
        The deconvolution results should be found at 
        adata.obsm[deconvolution_key]. 
        The default is 'cell_type_proportion'.
    spatial_info : str or list, optional
        The spatial information key in adata.obsm or columns in adata.obs. 
        The default is ['array_row', 'array_col'].
    shape : str
        The shape of the spots.
        The default is 'h'.
    dpi : int, optional
        Dots per inch. The default is 350.
    size : float list, optional
        Note there are two sizes of spots, corresponding to the large figures 
        and the small figures.
        The default is [7, 0.8].
    return_fig : bool, optional
        Where to return the figure. The default is False.
    save_path : system path, optional
        The path for the saving figure. The default is None.

    Returns
    -------
    fig : matplotlib axes
        Figure.

    """
    if 'ftu-' in ftu:
        ftu = ftu.replace('ftu-', '')
    gene_df = adata.var.loc[svg_list, :]

    ftu_total = [str(ind) for ind in range(1,
                                           1 + np.unique(gene_df['ftu']).size)]
    fig = plt.figure(dpi=dpi,
                     constrained_layout=True,
                     figsize=(12, 14))
    gene_df = gene_df.loc[svg_list, :]
    ftu_gene_list = gene_df.loc[gene_df['ftu'] == ftu,
                    :].index.tolist()

    # *************************************************
    # Spatial map plot scatter
    ax_SpaMap_scatter = plt.subplot2grid((12, 14), (1, 1), colspan=4,
                                         rowspan=4)
    ftu_spatial_map_scatter_ftu_id_card(adata,
                                        f"ftu_{ftu}",
                                        "#FF6879",
                                        radius=size[0],
                                        title=None,
                                        ax=ax_SpaMap_scatter,
                                        spatial_info=spatial_info,
                                        shape=shape,
                                        )
    # Spatial map plot start
    # Spatial map title
    ax_SpaMap_title = plt.subplot2grid((12, 14), (0, 1), colspan=4, rowspan=1)
    plt.title(f"Spatial map: ftu {ftu}", y=0, fontsize=20)
    ax_SpaMap_title.spines['top'].set_visible(False)
    ax_SpaMap_title.spines['right'].set_visible(False)
    ax_SpaMap_title.spines['left'].set_visible(False)
    ax_SpaMap_title.spines['bottom'].set_visible(False)
    ax_SpaMap_title.get_xaxis().set_visible(False)
    ax_SpaMap_title.get_yaxis().set_visible(False)

    # Spatial map plot end
    # *************************************************

    # *************************************************
    # Enhanced svgs plot start
    # Enhanced svgs title
    ax_enhsvgs_title = plt.subplot2grid((12, 14), (0, 6), colspan=4, rowspan=1)

    plt.title("Corresponding svgs", y=0, fontsize=20)
    ax_enhsvgs_title.spines['top'].set_visible(False)
    ax_enhsvgs_title.spines['right'].set_visible(False)
    ax_enhsvgs_title.spines['left'].set_visible(False)
    ax_enhsvgs_title.spines['bottom'].set_visible(False)
    ax_enhsvgs_title.get_xaxis().set_visible(False)
    ax_enhsvgs_title.get_yaxis().set_visible(False)
    # Enhanced svgs plot scatters
    for index in range(min(8, len(ftu_gene_list))):
        ax = plt.subplot2grid((12, 14), (1 + 2 * (index % 4),
                                         6 + 2 * (index // 4)),
                              colspan=2, rowspan=2)

        scatter_svgs_distri_ftu_id_card(adata,
                                        gene=ftu_gene_list[index],
                                        spatial_info=spatial_info,
                                        size=size[1],
                                        shape=shape,
                                        ax=ax)

    # Enhanced svgs plot end
    # *************************************************

    # *************************************************
    # Overlapped ftus plot start
    ftu_overlapper_max = {
        "ftu": None,
        "ob_sum": None
    }
    ftu_overlapper_min = {
        "ftu": None,
        "ob_sum": None
    }
    for value in ftu_total:
        v1 = [int(i) for i in adata.obsm["ftu_binary"][f"ftu_{value}"].values]
        v2 = [int(i) for i in adata.obsm["ftu_binary"][f"ftu_{ftu}"].values]
        current_value = [i for i in list(map(lambda x: x[0] + x[1],
                                             zip(v2, v1))) if i == 2]
        current_value = len(current_value)
        if ftu_overlapper_max["ftu"] == None and value != ftu:
            ftu_overlapper_max["ftu"] = value
            ftu_overlapper_max["ob_sum"] = current_value

        elif ftu_overlapper_min["ftu"] == None and value != ftu:
            ftu_overlapper_min["ftu"] = value
            ftu_overlapper_min["ob_sum"] = current_value
        elif value != ftu:
            if ftu_overlapper_min["ob_sum"] > current_value:
                ftu_overlapper_min["ob_sum"] = current_value
                ftu_overlapper_min["ftu"] = value
            if ftu_overlapper_max["ob_sum"] < current_value:
                ftu_overlapper_max["ob_sum"] = current_value
                ftu_overlapper_max["ftu"] = value
    ax_Overlapped_scatter_1 = plt.subplot2grid((12, 14), (5, 1),
                                               colspan=2,
                                               rowspan=2)
    ftu_overlapped_scatter_ftu_id_card(adata,
                                       ftu,
                                       ftu_overlapper_max["ftu"],
                                       ftu_1_color="#FF6879",
                                       ftu_2_color="Green",
                                       overlapped_color="Yellow",
                                       title="Overlapped ftus",
                                       radius=size[1],
                                       spatial_info=spatial_info,
                                       marker=shape,
                                       ax=ax_Overlapped_scatter_1)
    ax_Overlapped_scatter_2 = plt.subplot2grid((12, 14), (7, 1),
                                               colspan=2, rowspan=2)
    ftu_overlapped_scatter_ftu_id_card(adata,
                                       ftu,
                                       ftu_overlapper_min["ftu"],
                                       ftu_1_color="#FF6879",
                                       ftu_2_color="Green",
                                       overlapped_color="Yellow",
                                       title="Overlapped ftus",
                                       radius=size[1],
                                       marker=shape,
                                       spatial_info=spatial_info,
                                       ax=ax_Overlapped_scatter_2)
    ax_Overlapped_scatter_1_ftu = plt.subplot2grid((12, 14), (5, 3),
                                                   colspan=2, rowspan=2)

    ftu_spatial_map_scatter_ftu_id_card(adata,
                                        f"ftu_{str(ftu_overlapper_max['ftu'])}",
                                        "#FF6879",
                                        f"ftu {ftu_overlapper_max['ftu']}",
                                        ax=ax_Overlapped_scatter_1_ftu,
                                        radius=size[1],
                                        shape=shape,
                                        spatial_info=spatial_info)
    ax_Overlapped_scatter_2_ftu = plt.subplot2grid((12, 14), (7, 3),
                                                   colspan=2, rowspan=2)

    ftu_spatial_map_scatter_ftu_id_card(adata,
                                        f"ftu_{str(ftu_overlapper_min['ftu'])}",
                                        "#FF6879",
                                        f"ftu {ftu_overlapper_min['ftu']}",
                                        ax=ax_Overlapped_scatter_2_ftu,
                                        radius=size[1],
                                        shape=shape,
                                        spatial_info=spatial_info)
    # Overlapped ftus plot end
    # *************************************************

    # *************************************************
    # Cell type proportion plot start
    # Cell type proportion title
    cell_type_proportion_plot = plt.subplot2grid((12, 14), (7, 6),
                                                 colspan=4, rowspan=2)

    # Cell type proportion plot
    if deconvolution_key is not None:
        cell2loc = adata.obsm[deconvolution_key]
        sum_cell2loc = [sum(cell2loc[i].values.tolist()) for i in \
                        cell2loc.columns.tolist()]
        # print(sum_cell2loc)
        cell_type_name_i = [sum_cell2loc.index(x) for x in \
                            sorted(sum_cell2loc, reverse=True)[:10]]
        cell_type_name = [cell2loc.columns.tolist()[i] for i \
                          in cell_type_name_i]

        cell_type_proportion_box_ftu_id_card(cell_type_name,
                                             cell2loc,
                                             ax=cell_type_proportion_plot,
                                             title=deconvolution_key)
    # Cell type proportion plot end
    # *************************************************

    # *************************************************
    # svg functional enrichment plot start
    import gseapy as gp
    enr = gp.enrichr(gene_list=ftu_gene_list,
                     gene_sets=['BioPlanet_2019',
                                'GO_Biological_Process_2021',
                                'ChEA_2016'],
                     organism='Human',
                     description='ftu',
                     outdir='ftup/enrichr_kegg',
                     no_plot=False,
                     cutoff=0.5
                     )

    GO_Biological_Processn_plot = plt.subplot2grid((12, 14), (10, 1),
                                                   colspan=4, rowspan=2)

    from gseapy.plot import barplot
    GO_Biological_Process_file = "./ftup/enrichr_kegg/GO_Biological_Process"
    barplot(enr.results[enr.results.Gene_set == 'GO_Biological_Process_2021'],
            top_term=10, ofname=GO_Biological_Process_file, )
    from PIL import Image
    img = Image.open(GO_Biological_Process_file + ".png", )
    GO_Biological_Processn_plot.imshow(img)
    plt.title("GO Biological Process", fontsize=20)
    GO_Biological_Processn_plot.spines['top'].set_visible(False)
    GO_Biological_Processn_plot.spines['right'].set_visible(False)
    GO_Biological_Processn_plot.spines['left'].set_visible(False)
    GO_Biological_Processn_plot.spines['bottom'].set_visible(False)
    GO_Biological_Processn_plot.get_xaxis().set_visible(False)
    GO_Biological_Processn_plot.get_yaxis().set_visible(False)
    # svg functional enrichment's plot end
    # *************************************************

    # *************************************************

    ftu_heatmap_signal_plot = plt.subplot2grid((12, 14), (11, 6),
                                               colspan=2, rowspan=1)

    ftu_heatmap_signal_ftu_id_card(adata,
                                   ftu=f"ftu_{ftu}",
                                   ax=ftu_heatmap_signal_plot)

    ftu_freq_signal_plot = plt.subplot2grid((12, 14), (10, 6),
                                            colspan=2, rowspan=1)

    ftu_freq_signal_ftu_id_card(adata,
                                ftu=f"ftu_{ftu}",
                                ax=ftu_freq_signal_plot,
                                title="FCs")
    logo_plot = plt.subplot2grid((12, 14), (10, 8), colspan=2, rowspan=2)

    # Logo_plot_file = "./SpaGFT_Logo_RGB.png"
    # from PIL import Image
    # img = Image.open(Logo_plot_file)
    # Logo_plot.imshow(img)
    logo_plot.spines['top'].set_visible(False)
    logo_plot.spines['right'].set_visible(False)
    logo_plot.spines['left'].set_visible(False)
    logo_plot.spines['bottom'].set_visible(False)
    logo_plot.get_xaxis().set_visible(False)
    logo_plot.get_yaxis().set_visible(False)
    plt.text(x=0, y=0.5, s="Powered \n      by  \n SpaGFT",
             fontsize=20, fontstyle="italic",
             fontweight="bold")
    plt.title("ftu " + ftu + " ID Card \n# of svgs: " + \
              str(len(ftu_gene_list)), y=0)
    # # *************************************************
    # plt.tight_layout()
    if save_path is not None:
        plt.savefig(f"{save_path}/ftu_{ftu}.png")
    if os.path.exists("./ftup/enrichr_kegg"):
        os.system("rm -r ./ftup/enrichr_kegg")
    if return_fig:
        return fig
    plt.show()
    plt.close()
