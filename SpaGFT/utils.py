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
        adata.X is the normalized count matrix. Besides, the spatial coordinates could be found in adata.obs
        or adata.obsm.
    num_neighbors: int, optional
        The number of neighbors for each node/spot/pixel when construct graph.
        The default if 6.
    spatial_key=None : list | string
        Get the coordinate information by adata.obsm[spaital_key] or adata.var[spatial_key].
        The default is ['array_row', 'array_col'].
    normalization : bool, optional
        Whether you need to normalize laplacian matrix. The default is False.

    Raises
    ------
    KeyError
        The coordinates should be found at adata.obs[spatial_names] or adata.obsm[spatial_key]

    Returns
    -------
    lap_mtx : csr_matrix
        The Laplacian matrix or normalized Laplacian matrix.

    """
    if spatial_key in adata.obsm_keys():
        adj_mtx = kneighbors_graph(adata.obsm[spatial_key],
                                   n_neighbors=num_neighbors)
    elif set(spatial_key) <= set(adata.obs_keys()):
        adj_mtx = kneighbors_graph(adata.obs[spatial_key],
                                   n_neighbors=num_neighbors)
    else:
        raise KeyError("%s is not available in adata.obsm_keys" % spatial_key + " or adata.obs_keys")

    adj_mtx = nx.adjacency_matrix(nx.Graph(adj_mtx))
    # Degree matrix
    deg_mtx = adj_mtx.sum(axis=1)
    deg_mtx = create_degree_mtx(deg_mtx)
    # Laplacian matrix
    # Whether you need to normalize Laplacian matrix
    if not normalization:
        lap_mtx = deg_mtx - adj_mtx
    else:
        deg_mtx = np.array(adj_mtx.sum(axis=1)) ** (-0.5)
        deg_mtx = create_degree_mtx(deg_mtx)
        lap_mtx = ss.identity(deg_mtx.shape[0]) - deg_mtx @ adj_mtx @ deg_mtx

    return lap_mtx


def create_adjacent_mtx(coor_df,
                        spatial_names=['array_row', 'array_col'],
                        num_neighbors=4):
    # Transform coordinate dataframe to coordinate array
    coor_array = coor_df.loc[:, spatial_names].values
    coor_array.astype(np.float32)
    edge_list = []
    num_neighbors += 1
    for i in range(coor_array.shape[0]):
        point = coor_array[i, :]
        distances = np.sum(np.asarray((point - coor_array) ** 2), axis=1)
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
    # Normalization
    frequency_array = preprocessing.StandardScaler().fit_transform(frequency_array)

    # Dimension reduction
    if reduction and frequency_array.shape[1] > reduction:
        pca = PCA(n_components=reduction)
        frequency_array = pca.fit_transform(frequency_array)

    # Clustering
    kmeans_model = KMeans(n_clusters=n_clusters).fit(frequency_array)

    return kmeans_model.labels_


def window_side_bin(adata,
                    shape=(20, 20),
                    spatial_names=['array_row', 'array_col'],
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
                count_mtx.loc["pseudo_" + str(i) + "_" + str(j), :] = adata[tmp_index, :].X.sum(axis=0)
                final_coor_df.loc["pseudo_" + str(i) + "_" + str(j), :] = [i, j]
            else:
                count_mtx.loc["pseudo_" + str(i) + "_" + str(j), :] = 0
                final_coor_df.loc["pseudo_" + str(i) + "_" + str(j), :] = [i, j]

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


def kneed_select_values(value_list, S=3, increasing=True):
    from kneed import KneeLocator
    x_list = list(range(1, 1 + len(value_list)))
    y_list = value_list.copy()
    if increasing:
        magic = KneeLocator(x=x_list,
                            y=y_list,
                            S=S)
    else:
        y_list = y_list[::-1].copy()
        magic = KneeLocator(x=x_list,
                            y=y_list,
                            direction='decreasing',
                            S=S,
                            curve='convex')
    return magic.elbow


def correct_pvalues_for_multiple_testing(pvalues,
                                         correction_type="Benjamini-Hochberg"):
    """
    Correct p-values to obtain the adjusted p-values

    Parameters
    ----------
    pvalues : list | 1-D array
        The original p values. It should be a list.
    correction_type : str, optional
        Method used to correct p values. The default is "Benjamini-Hochberg".

    Returns
    -------
    new_pvalues : array
        Corrected p values.

    """
    from numpy import array, empty
    pvalues = array(pvalues)
    n = int(pvalues.shape[0])
    new_pvalues = empty(n)
    if correction_type == "Bonferroni":
        new_pvalues = n * pvalues
    elif correction_type == "Bonferroni-Holm":
        values = [(pvalue, i) for i, pvalue in enumerate(pvalues)]
        values.sort()
        for rank, vals in enumerate(values):
            pvalue, i = vals
            new_pvalues[i] = (n - rank) * pvalue
    elif correction_type == "Benjamini-Hochberg":
        values = [(pvalue, i) for i, pvalue in enumerate(pvalues)]
        values.sort()
        values.reverse()
        new_values = []
        for i, vals in enumerate(values):
            rank = n - i
            pvalue, index = vals
            new_values.append((n / rank) * pvalue)
        for i in range(0, int(n) - 1):
            if new_values[i] < new_values[i + 1]:
                new_values[i + 1] = new_values[i]
        for i, vals in enumerate(values):
            pvalue, index = vals
            new_pvalues[index] = new_values[i]
    return new_pvalues


def permutation_signal(signal_array, num_permutation=1000):
    """
    Permutate gene signals in spatial domain randomly.

    Parameters
    ----------
    signal_array : list | array
        A one-dimensional array indicate gene expression on all spots.
    num_permutation : int, optional
        The number of permutation. The default is 1000.

    Returns
    -------
    total_signals : array
        The permuted gene expression signals.

    """
    signal_array = np.array(signal_array)
    total_signals = signal_array * np.ones((num_permutation,
                                            len(signal_array)))

    for i in range(num_permutation):
        total_signals[i, :] = np.random.permutation(total_signals[i, :])

    return total_signals


def significant_test_permutation(exp_mtx,
                                 gene_score,
                                 eigvals,
                                 eigvecs_T,
                                 num_permutaion=1000,
                                 num_pool=200,
                                 spec_norm='l1'):
    """
    To calculate p values for genes, permutate gene expression data and 
    calculate p values.

    Parameters
    ----------
    exp_mtx : 2D-array
        The count matrix of gene expressions. (spots * genes)
    gene_score : 1D-array
        The calculated gene scores. 
    eigvals : array
        The eigenvalues of Laplacian matrix.
    eigvecs_T : array
        The eigenvectors of Laplacian matrix.
    num_permutation : int, optional
        The number of permutations. The default is 1000.
    num_pool : int, optional
        The cores used for multiprocess calculation to accelerate speed. The default is 200.
    spec_norm : str, optional
        The method to normalize graph signals in spectral domain. The default  is 'l1'.

    Returns
    -------
    array
        The calculated p values.

    """
    from multiprocessing.dummy import Pool as ThreadPool
    from scipy.stats import mannwhitneyu

    def _test_by_permutaion(gene_index):
        (gene_index)
        graph_signal = exp_mtx[gene_index, :]
        total_signals = permutation_signal(signal_array=graph_signal, num_permutation=num_permutaion)
        frequency_array = np.matmul(eigvecs_T, total_signals.transpose())
        frequency_array = np.abs(frequency_array)
        if spec_norm != None:
            frequency_array = preprocessing.normalize(frequency_array,
                                                      norm=spec_norm,
                                                      axis=0)
        score_list = np.matmul(2 ** (-1 * eigvals), frequency_array)
        score_list = score_list / score_max
        pval = mannwhitneyu(score_list, gene_score[gene_index], alternative='less').pvalue
        return pval

    score_max = np.matmul(2 ** (-2 * eigvals), (1 / len(eigvals)) * \
                          np.ones(len(eigvals)))
    gene_index_list = list(range(exp_mtx.shape[0]))
    pool = ThreadPool(num_pool)
    res = pool.map(_test_by_permutaion, gene_index_list)

    return res


def test_significant_freq(freq_array,
                          cutoff,
                          num_pool=200):
    """
    Significance test by comparing the intensities in low frequency FMs and in high frequency FMs.

    Parameters
    ----------
    freq_array : array
        The graph signals of genes in frequency domain. 
    cutoff : int
        Watershed between low frequency signals and high frequency signals.
    num_pool : int, optional
        The cores used for multiprocess calculation to accelerate speed. The
        default is 200.

    Returns
    -------
    array
        The calculated p values.

    """
    from scipy.stats import wilcoxon, mannwhitneyu, ranksums, combine_pvalues
    from multiprocessing.dummy import Pool as ThreadPool

    def _test_by_feq(gene_index):
        freq_signal = freq_array[gene_index, :]
        freq_1 = freq_signal[:cutoff]
        freq_1 = freq_1[freq_1 > 0]
        freq_2 = freq_signal[cutoff:]
        freq_2 = freq_2[freq_2 > 0]
        if freq_1.size <= 80 or freq_2.size <= 80:
            freq_1 = np.concatenate((freq_1, freq_1, freq_1, freq_1))
            freq_2 = np.concatenate((freq_2, freq_2, freq_2, freq_2))
        if freq_1.size <= 120 or freq_2.size <= 120:
            freq_1 = np.concatenate((freq_1, freq_1, freq_1))
            freq_2 = np.concatenate((freq_2, freq_2, freq_2))
        if freq_1.size <= 160 or freq_2.size <= 160:
            freq_1 = np.concatenate((freq_1, freq_1))
            freq_2 = np.concatenate((freq_2, freq_2))
        pval = ranksums(freq_1, freq_2, alternative='greater').pvalue
        return pval

    gene_index_list = list(range(freq_array.shape[0]))
    pool = ThreadPool(num_pool)
    res = pool.map(_test_by_feq, gene_index_list)

    return res


def my_eigsh(args_tuple):
    """
    The function is used to multi-process calculate using Pool.

    Parameters
    ----------
    args_tuple : tupple
        The args_tupple contains three elements, that are, Lplacian matrix, k 
        and which.

    Returns
    -------
    (eigvals, eigvecs)

    """
    lap_mtx = args_tuple[0]
    k = args_tuple[1]
    which = args_tuple[2]
    eigvals, eigvecs = ss.linalg.eigsh(lap_mtx.astype(float),
                                       k=k,
                                       which=which)
    return ((eigvals, eigvecs))


def get_cos_similar(v1: list, v2: list):
    v3 = np.array(v1) + np.array(v2)
    return v3[v3 >= 2].size


def get_overlap_cs_core(cluster_collection):
    over_loop_cs_score = 0
    over_loop_cs_score_num = 0
    spot_shape = cluster_collection.shape[1]
    for index, _ in enumerate(cluster_collection):
        v1 = cluster_collection[index]
        if index + 1 <= cluster_collection.shape[0]:
            for value in cluster_collection[index + 1:]:
                over_loop_cs_score += get_cos_similar(value, v1)
                over_loop_cs_score_num += 1
    return (over_loop_cs_score / over_loop_cs_score_num / spot_shape) \
        if over_loop_cs_score_num != 0 else 1
