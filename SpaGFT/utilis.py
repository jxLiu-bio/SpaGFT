import numpy as np
import pandas as pd
import scanpy as sc
from scipy.spatial.distance import pdist, squareform
import itertools
import scipy.sparse as ss
from sklearn import preprocessing
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
from sklearn.neighbors import kneighbors_graph
import networkx as nx


def get_laplacian_mtx(adata, 
                      num_neighbors=6,
                      spatial_key=['array_row', 'array_col'],
                      normalization=False):
    """
    Obtain the Laplacian matrix or normalized laplacian matrix.

    Parameters
    ----------
    adata : AnnData
        adata.X is the normalized count matrix. Besides, the spatial coordinat-
        es could be found in adata.obs or adata.obsm.
    num_neighbors: int, optional
        The number of neighbors for each node/spot/pixel when contrcut graph. 
        The defalut if 6.
    spatial_key=None : list | string
        Get the coordinate information by adata.obsm[spaital_key] or 
        adata.var[spatial_key]. The default is ['array_row', 'array_col'].
    normalization : bool, optional
        Whether need to normalize laplacian matrix. The default is False.

    Raises
    ------
    KeyError
        The coordinates should be found at adata.obs[spatial_names] or 
        adata.obsm[spatial_key]

    Returns
    -------
    lap_mtx : csr_matrix
        The laplcaian matrix or mormalized laplcian matrix.

    """
    if spatial_key in adata.obsm_keys():
        adj_mtx = kneighbors_graph(adata.obsm[spatial_key], 
                                    n_neighbors=num_neighbors)
    elif set(spatial_key) <= set(adata.obs_keys()):
        adj_mtx = kneighbors_graph(adata.obs[spatial_key], 
                                    n_neighbors=num_neighbors)
    else:
        raise KeyError("%s is not avaliable in adata.obsm_keys" % \
                       spatial_key +  " or adata.obs_keys")

    adj_mtx = nx.adjacency_matrix(nx.Graph(adj_mtx))
    # Degree matrix
    deg_mtx = adj_mtx.sum(axis=1)
    deg_mtx = create_degree_mtx(deg_mtx)
    # Laplacian matrix
    # Whether need to normalize laplcian matrix
    if not normalization:
        lap_mtx = deg_mtx - adj_mtx
    else:
        deg_mtx = np.array(adj_mtx.sum(axis=1)) ** (-0.5)
        deg_mtx = create_degree_mtx(deg_mtx)
        lap_mtx = ss.identity(deg_mtx.shape[0]) - deg_mtx @ adj_mtx @ deg_mtx
    
    return lap_mtx

def find_HVGs(adata, norm_method=None, num_genes=2000):
    """
    Find spatialy variable genes.

    Parameters
    ----------
    adata : AnnData
        The object to save sequencing data.
    norm_method : str | None, optional
        The method used to normalized adata. The default is None.
    num_genes : None | int, optional
        The number of highly variable genes. The default is 2000.

    Returns
    -------
    HVG_list : list
        Highly variable genes.

    """
    # Normalization
    if norm_method == "CPM":
        sc.pp.normalize_total(adata, target_sum=1e5)
        # log-transform
        sc.pp.log1p(adata)
    else:
        pass
    # Find high variable genes using sc.pp.highly_variable_genes() with default 
    # parameters
    sc.pp.highly_variable_genes(adata, n_top_genes=num_genes)
    HVG_list = adata.var.index[adata.var.highly_variable]
    HVG_list = list(HVG_list)
    
    return HVG_list

def create_adjacent_mtx(coor_df, spatial_names = ['array_row', 'array_col'],
                        num_neighbors=4):
    # Transform coordinate dataframe to coordinate array
    coor_array = coor_df.loc[:, spatial_names].values
    coor_array.astype(np.float32) 
    edge_list = []
    num_neighbors += 1
    for i in range(coor_array.shape[0]):
        point = coor_array[i, :]
        distances = np.sum(np.asarray((point - coor_array)**2), axis=1)
        distances = pd.DataFrame(distances, 
                                 index=range(coor_array.shape[0]),
                                 columns=["distance"])
        distances = distances.sort_values(by='distance', ascending=True)
        neighbors = distances[1:num_neighbors].index.tolist()
        edge_list.extend((i, j) for j in neighbors)
        edge_list.extend((j, i) for j in neighbors)
    # Remove duplicates
    edge_list = set(edge_list)
    edge_list = list(edge_list)
    row_index = []
    col_index = []
    row_index.extend(j[0] for j in edge_list)
    col_index.extend(j[1] for j in edge_list)
    
    sparse_mtx = ss.coo_matrix((np.ones_like(row_index), (row_index, col_index)),
                               shape=(coor_array.shape[0], coor_array.shape[0]))
    
    return sparse_mtx

def create_degree_mtx(diag):
    diag = np.array(diag)
    diag = diag.flatten()
    row_index = list(range(diag.size))
    col_index = row_index
    sparse_mtx = ss.coo_matrix((diag, (row_index, col_index)),
                               shape=(diag.size, diag.size))
    
    return sparse_mtx

def gene_clustering_kMeans(frequency_array, n_clusters, reduction=50):
    # Normlization
    frequency_array = preprocessing.StandardScaler().fit_transform(frequency_array)
    
    # Dimension reduction
    if reduction and frequency_array.shape[1] > reduction:
        pca = PCA(n_components=reduction)
        frequency_array = pca.fit_transform(frequency_array)
            
    # Clustering
    kmeans_model = KMeans(n_clusters=n_clusters).fit(frequency_array)
    
    return kmeans_model.labels_

def window_side_bin(adata, shape=(20, 20),
                    spatial_names = ['array_row', 'array_col'],
                    sparse=True):
    # Extract border
    max_x = adata.obs[spatial_names[1]].max()   
    min_x = adata.obs[spatial_names[1]].min()
    max_y = adata.obs[spatial_names[0]].max()
    min_y = adata.obs[spatial_names[0]].min()
    
    # Calculate bin-size
    bin_x = (max_x - min_x) / (shape[0] - 1)
    bin_y = (max_y - min_y) / (shape[1] - 1)
    

    # Create a new dataframe to store new coordinates
    new_coor_df = pd.DataFrame(0, index=adata.obs.index, columns=spatial_names)
    for i in range(adata.shape[0]):
        coor_x = adata.obs.iloc[i][spatial_names][1]
        coor_y = adata.obs.iloc[i][spatial_names][0]
        coor_x = int(np.floor((coor_x - min_x) / bin_x))
        coor_y = int(np.floor((coor_y - min_y) / bin_y))
        new_coor_df.iloc[i, :] = [coor_x, coor_y]
    
    # Merge the spots within the bins
    count_mtx = pd.DataFrame(columns=adata.var_names)
    final_coor_df = pd.DataFrame(columns=spatial_names)
    for i in range(shape[0]):
        for j in range(shape[1]):
            tmp_index = new_coor_df[new_coor_df.iloc[:, 0] == i]
            tmp_index = tmp_index[tmp_index.iloc[:, 1] == j]
            if tmp_index.shape[0] > 0:
                tmp_index = tmp_index.index
                count_mtx.loc["pseudo_" + str(i) + "_" + str(j), 
                              :] = adata[tmp_index, :].X.sum(axis=0)
                final_coor_df.loc["pseudo_" + str(i) + "_" + str(j),
                                  :] = [i, j]
            else:
                count_mtx.loc["pseudo_" + str(i) + "_" + str(j), 
                              :] = 0
                final_coor_df.loc["pseudo_" + str(i) + "_" + str(j),
                                  :] = [i, j]
    
    # Transform it to anndata
    from anndata import AnnData
    if sparse == True:
        from scipy import sparse
        new_adata = AnnData(sparse.coo_matrix(count_mtx))
        new_adata.obs_names = count_mtx.index.tolist()
        new_adata.var_names = count_mtx.columns.tolist()
        new_adata.obs = final_coor_df
    else:
        new_adata = AnnData(count_mtx)  
        new_adata.obs = final_coor_df
        
    return new_adata

def select_svg_normal(gene_score, num_sigma=1):
    mu = np.mean(gene_score['gft_score'])
    sigma = np.std(gene_score['gft_score'])
    gene_score['spatially_variable'] = 0
    gene_score.loc[gene_score['gft_score'] > mu + num_sigma * sigma, 
                   'spatially_variable'] = 1
    
    return gene_score

def select_svg_kmean(gene_score):
    from sklearn.cluster import KMeans
    kmeans = KMeans(n_clusters=2)
    X = gene_score["smooth_score"].tolist()
    X = np.array(X).reshape(-1, 1)
    y_pred = kmeans.fit_predict(X)       
    gene_score['spatially_variable'] = y_pred
    
    return gene_score

def umap_spectral_domain(frequency_array, gene_score, n_dim=2):
    adata_gene = sc.AnnData(frequency_array).T
    sc.pp.pca(adata_gene)
    sc.pp.neighbors(adata_gene)
    sc.tl.umap(adata_gene)
    gene_score = select_svg_normal(gene_score, num_sigma=1)
    gene_score = gene_score.reindex(adata_gene.obs.index)
    adata_gene.obs = gene_score
    sc.pl.umap(adata_gene, color='spatially_variable')
    
def tsne_spectral_domain(frequency_array, gene_score, n_dims=2):
    adata_gene = sc.AnnData(frequency_array).T
    sc.pp.pca(adata_gene)
    sc.pp.neighbors(adata_gene)
    sc.tl.tsne(adata_gene)
    gene_score = select_svg_normal(gene_score, num_sigma=1)
    gene_score = gene_score.reindex(adata_gene.obs.index)
    adata_gene.obs = gene_score
    sc.pl.tsne(adata_gene, color='spatially_variable')

def fms_spectral_domain(frequency_array, gene_score, n_dims=2):
    adata_gene = sc.AnnData(frequency_array).T
    adata_gene.obsm['X_pca'] = frequency_array.transpose()
    gene_score = select_svg_normal(gene_score, num_sigma=3)
    gene_score = gene_score.reindex(adata_gene.obs.index)
    adata_gene.obs = gene_score
    sc.pl.pca(adata_gene, color='spatially_variable')

def pca_spatial_domain(adata, gene_score, n_dims=2):
    adata_gene = adata.copy()
    adata_gene = adata_gene.T
    gene_score = select_svg_normal(gene_score, num_sigma=1)
    adata_gene.obs['spatially_variable'] = 0
    svg_index = (gene_score[gene_score['spatially_variable'] == 1]).index
    adata_gene.obs.loc[svg_index, 'spatially_variable'] = 1
    sc.pp.pca(adata_gene)
    sc.pl.pca(adata_gene, color='spatially_variable')

def umap_spatial_domain(adata, gene_score, n_dim=2):
    adata_gene = adata.copy()
    adata_gene = adata_gene.T
    gene_score = select_svg_normal(gene_score, num_sigma=1)
    adata_gene.obs['spatially_variable'] = 0
    svg_index = (gene_score[gene_score['spatially_variable'] == 1]).index
    adata_gene.obs.loc[svg_index, 'spatially_variable'] = 1
    sc.pp.pca(adata_gene)
    sc.pp.neighbors(adata_gene)
    sc.tl.umap(adata_gene)
    sc.pl.umap(adata_gene, color='spatially_variable')

def tsne_spatial_domain(adata, gene_score, n_dim=2):
    adata_gene = adata.copy()
    adata_gene = adata_gene.T
    gene_score = select_svg_normal(gene_score, num_sigma=1)
    adata_gene.obs['spatially_variable'] = 0
    svg_index = (gene_score[gene_score['spatially_variable'] == 1]).index
    adata_gene.obs.loc[svg_index, 'spatially_variable'] = 1
    sc.pp.pca(adata_gene)
    sc.pp.neighbors(adata_gene)
    sc.tl.tsne(adata_gene)
    sc.pl.tsne(adata_gene, color='spatially_variable')
    
    
def cal_mean_expression(adata, gene_list):
    tmp_adata = adata[:, gene_list].copy()
    if 'log1p' not in adata.uns_keys():
        tmp_adata = sc.pp.log1p(tmp_adata)
    mean_vector = tmp_adata.X.mean(axis=1)
    mean_vector = np.array(mean_vector).ravel()
    
    return mean_vector





    
                
                
        
        
    
            
        
        
        
        
        
        
    
    

