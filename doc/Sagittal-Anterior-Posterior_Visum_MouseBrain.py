import numpy as np
import pandas as pd
import scanpy as sc
import SpaGFT as spg
import os
import matplotlib.pyplot as plt
from plotnine import *
import gseapy as gp


sc.settings.set_figure_params(dpi=100, frameon=False, figsize=(7, 7),
                              facecolor='white')
# Load data and preprocessing
data_dir = "/home/frank/Desktop/SpaGFT_git/SpaGFT_codes/datasets"
sample_list = ["Sagittal-Anterior-1_Visium_MouseBrain",
               "Sagittal-Anterior-2_Visium_MouseBrain",
               "Sagittal-Posterior-1_Visium_MouseBrain",
               "Sagittal-Posterior-2_Visium_MouseBrain"]
for sample in sample_list:
    adata = sc.read_visium(os.path.join(data_dir, sample))
    adata.var_names_make_unique()
    adata.raw = adata
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
                                   filter_peaks=True,
                                   S=5)
    svg_list = gene_df[gene_df.cutoff_gft_score]\
        [gene_df.qvalue < 0.05].index.tolist()
    gene_df = gene_df.loc[svg_list, :]
    # Visualize some SVGs and their spatial domain
    for gene in svg_list[:2]:
        sc.pl.spatial(adata, color=gene, cmap='magma')
        spg.plot.svg_freq_signal(adata, gene)
    # Visualize all genes' frequency signals by UMAP
    spg.plot.gene_signal_umap(adata, svg_list, s=30)
    
    
    # *************************************
    # Detect tissue modules
    # *************************************
    spg.gft.find_tissue_module(adata,
                               svg_list = svg_list,
                               ratio_fms=2,
                               resolution=1,
                               n_neighbors=15,
                               sub_resolution=0.5,
                               sub_n_neighbors=8,
                               quantile=0.75)
    gene_df = adata.var.loc[svg_list, :]
    # UMAP figure of clustering result
    plot_df = pd.concat((adata.uns['gft_umap_tm'], gene_df.tissue_module), axis=1)
    plot_df.tissue_module = pd.Categorical(plot_df.tissue_module,
                                    categories=np.unique(plot_df.tissue_module))
    base_plot = (ggplot(plot_df, aes('UMAP_1', 'UMAP_2', fill='tissue_module'))
                  + geom_point(size=3.5, stroke=0)
                  + scale_fill_hue(s=0.99, l=0.65, h=0.0417, color_space='husl')
                  + theme_classic())
    print(base_plot)
    # Visualize some genes determining a tissue module
    for gene in gene_df[gene_df.tissue_module == '7'].index[:3]:
        sc.pl.spatial(adata, color=gene)
    # Tissue modules
    tm_df = adata[:, svg_list].var
    # Visualize frequency signals of some Tissue Modules or sub-TM
    all_tms = adata.uns['freq_signal_tm'].index[:3]
    for tm in all_tms:
        signal = adata.uns['freq_signal_tm'].loc[tm, :]
        ax = spg.plot.tm_freq_signal(adata, tm, return_fig=True)
    all_subTMs = adata.uns['freq_signal_subTM'].index[:3]
    for subTM in all_subTMs:
        signal = adata.uns['freq_signal_subTM'].loc[subTM, :]
        ax = spg.plot.subTm_freq_signal(adata, subTM, return_fig=True)
    # Visualize on tissue module, including pseudo_expression
    adata.obs['tm_1'] = adata.obsm['tm_pseudo_expression']['tm_1']
    sc.pl.spatial(adata, color='tm_1')
    adata.obs['tm_1_binary'] = adata.obsm['tm_binary']['tm_1']
    sc.pl.spatial(adata, color='tm_1_binary')
    adata.uns['tm_1_binary_colors'] = ['#CCCCCC', '#CA1C1C']
    sc.pl.spatial(adata)
    sc.pl.spatial(adata, color='tm_1_binary')
    # Biological pathway enrichment analysis for detected tissue modules
    tm_gene_list = tm_df[tm_df.tissue_module == '1'].index.tolist()
    enr = gp.enrichr(gene_list=tm_gene_list,
                      gene_sets=['BioPlanet_2019','GO_Biological_Process_2021',
                                'ChEA_2016'],
                      organism='Human', 
                      description='Tissue_module',
                      outdir='test/enrichr_kegg',
                      no_plot=False,
                      cutoff=0.5 # test dataset, use lower value from range(0,1)
                    )
    enr_results = enr.results
    from gseapy.plot import barplot, dotplot
    barplot(enr.results[enr.results.Gene_set=='GO_Biological_Process_2021'],
            title='GO_Biological_Process_2021')
    barplot(enr.results[enr.results.Gene_set=='BioPlanet_2019'],
            title='BioPlanet_2019')
    
    
    # **************************
    # Enhancement
    # **************************
    new_adata = adata.copy()
    spg.low_pass_enhancement(new_adata,
                            ratio_low_freq=15,
                            inplace=True)
    # Before
    for gene in svg_list[400:402]:
        sc.pl.spatial(adata, color=gene, cmap='magma', use_raw=False)
    # After
    for gene in svg_list[400:402]:
        sc.pl.spatial(new_adata, color=gene, cmap='magma', use_raw=False)
    
    
    
    
