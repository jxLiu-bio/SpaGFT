import warnings
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as ss
from sklearn import preprocessing
from sklearn.cluster import KMeans

from SpaGFT.utils import get_laplacian_mtx, kneed_select_values
from SpaGFT.utils import test_significant_freq, get_overlap_cs_core

warnings.filterwarnings("ignore")


def low_pass_enhancement(adata,
                         ratio_low_freq='infer',
                         ratio_neighbors='infer',
                         c=0.001,
                         spatial_info=['array_row', 'array_col'],
                         normalize_lap=False,
                         inplace=False):
    """
    Implement gene expressions with low-pass filter. After this step, the spatially variables genes will be smoother
    than the previous. The function can also be treated as denoising function. Note that the denoising results
    is related to spatial graph topology so that only the results of spatially variable genes are reasonable.

    Parameters
    ----------
    adata : AnnData
        adata.X is the normalized count matrix. Besides, the spatial coordinates of all spots should be found in
        adata.obs or adata.obsm.
    ratio_low_freq : float | "infer", optional
        The ratio_low_freq will be used to determine the number of the FMs with low frequencies. Indeed,
        the ratio_low_freq * sqrt(number of spots) low frequency FMs will be calculated. The default is 'infer'.
    ratio_neighbors: float | 'infer', optional
        The ratio_neighbors will be used to determine the number of neighbors
        when construct the KNN graph by spatial coordinates. Indeed, ratio_neighbors * sqrt(number of spots) / 2
        indicates the K. If 'infer', the parameter will be set to 1.0. The default is 'infer'.
    c: float, optional
        c balances the smoothness and difference with previous expression. A high can achieve better smoothness.
        c should be set to [0, 0.1].The default is 0.001.
    spatial_info : list | tuple | string, optional
        The column names of spatial coordinates in adata.obs_names or key in adata.obsm_keys() to obtain spatial
        information. The default is ['array_row', 'array_col'].
    normalize_lap : bool. optional
        Whether you need to normalize the Laplacian matrix. The default is False.
    inplace: bool, optional
        Whether you need to replace adata.X with the enhanced expression matrix. The default is False.
        

    Returns
    -------
    adata: anndata

    """
    import scipy.sparse as ss
    if ratio_low_freq == 'infer':
        if adata.shape[0] <= 800:
            num_low_frequency = min(20 * int(np.ceil(np.sqrt(adata.shape[0]))),
                                    adata.shape[0])
        elif adata.shape[0] <= 5000:
            num_low_frequency = 15 * int(np.ceil(np.sqrt(adata.shape[0])))
        elif adata.shape[0] <= 10000:
            num_low_frequency = 10 * int(np.ceil(np.sqrt(adata.shape[0])))
        else:
            num_low_frequency = 5 * int(np.ceil(np.sqrt(adata.shape[0])))
    else:
        num_low_frequency = int(np.ceil(np.sqrt(adata.shape[0]) * ratio_low_freq))

    if ratio_neighbors == 'infer':
        if adata.shape[0] <= 500:
            num_neighbors = 4
        else:
            num_neighbors = int(np.ceil(np.sqrt(adata.shape[0]) / 2))
    else:
        num_neighbors = int(np.ceil(np.sqrt(adata.shape[0]) / 2 * ratio_neighbors))

    adata.var_names_make_unique()
    sc.pp.filter_genes(adata, min_cells=1)

    # Get Laplacian matrix according to coordinates 
    lap_mtx = get_laplacian_mtx(adata,
                                num_neighbors=num_neighbors,
                                spatial_key=spatial_info,
                                normalization=normalize_lap)

    # Fourier mode of low frequency
    num_low_frequency = min(num_low_frequency, adata.shape[0])
    eigvals, eigvecs = ss.linalg.eigsh(lap_mtx.astype(float),
                                       k=num_low_frequency,
                                       which='SM')

    # *********************** Graph Fourier Transform **************************
    # Calculate GFT
    eigvecs_t = eigvecs.transpose()
    if not ss.issparse(adata.X):
        exp_mtx = adata.X
    else:
        exp_mtx = adata.X.toarray()
    frequency_array = np.matmul(eigvecs_t, exp_mtx)
    # low-pass filter
    filter_list = [1 / (1 + c * eigv) for eigv in eigvals]
    # filter_list = [np.exp(-c * eigv) for eigv in eigvals]
    filter_array = np.matmul(np.diag(filter_list), frequency_array)
    filter_array = np.matmul(eigvecs, filter_array)
    filter_array[filter_array < 0] = 0

    # whether you need to replace original count matrix
    if inplace and not ss.issparse(adata.X):
        adata.X = filter_array
    elif inplace:
        import scipy.sparse as ss
        adata.X = ss.csr.csr_matrix(filter_array)

    return adata


def determine_frequency_ratio(adata,
                              low_end=5,
                              high_end=5,
                              ratio_neighbors='infer',
                              spatial_info=['array_row', 'array_col'],
                              normalize_lap=False):
    '''
    This function can choose the number of FMs automatically based on
    kneedle algorithm.

    Parameters
    ----------
    adata : AnnData
        adata.X is the normalized count matrix. Besides, the spatial coordinates of all spots should be found in
        adata.obs or adata.obsm.
    low_end : float, optional
        The range of low-frequency FMs. The default is 5. Note that the real cutoff is low_end * sqrt(n_spots).
    high_end : TYPE, optional
        The range of high-frequency FMs. The default is 5.Note that the real cutoff is low_end * sqrt(n_spots).
    ratio_neighbors : float, optional
        The ratio_neighbors will be used to determine the number of neighbors when construct the KNN graph by spatial
        coordinates. Indeed, ratio_neighbors * sqrt(number of spots) / 2 indicates the K. If 'infer', the parameter
        will be set to 1.0. The default is 'infer'.
    spatial_info : list | tuple | string, optional
        The column names of spatial coordinates in adata.obs_names or key
        in adata.obsm_keys() to obtain spatial information. The default
        is ['array_row', 'array_col'].
    normalize_lap : bool, optional
        Whether you need to normalize the Laplacian matrix. The default is False.

    Returns
    -------
    low_cutoff : float
        The low_cutoff * sqrt(the number of spots) low-frequency FMs are 
        recommended in detecting svg.
    high_cutoff : float
        The high_cutoff * sqrt(the number of spots) low-frequency FMs are 
        recommended in detecting svg.

    '''
    # Determine the number of neighbors
    if ratio_neighbors == 'infer':
        num_neighbors = int(np.ceil(np.sqrt(adata.shape[0]) / 2))
    else:
        num_neighbors = int(np.ceil(np.sqrt(adata.shape[0]) / 2 \
                                    * ratio_neighbors))
    if adata.shape[0] <= 500:
        num_neighbors = 4
    if adata.shape[0] > 15000 and low_end >= 3:
        low_end = 3
    if adata.shape[0] > 15000 and high_end >= 3:
        high_end = 3
    # Ensure gene index uniquely and all gene had expression  
    adata.var_names_make_unique()
    sc.pp.filter_genes(adata, min_cells=1)

    # *************** Construct graph and corresponding matrixs ***************
    lap_mtx = get_laplacian_mtx(adata,
                                num_neighbors=num_neighbors,
                                spatial_key=spatial_info,
                                normalization=normalize_lap)
    print("Obatain the Laplacian matrix")

    # Next, calculate the eigenvalues and eigenvectors of the Laplace matrix
    # Fourier bases of low frequency
    eigvals_s, eigvecs_s = ss.linalg.eigsh(lap_mtx.astype(float),
                                           k=int(np.ceil(low_end * np.sqrt(adata.shape[0]))),
                                           which='SM')
    low_cutoff = np.ceil(kneed_select_values(eigvals_s) / np.sqrt(adata.shape[0]) * 1000) / 1000
    if low_cutoff >= low_end:
        low_cutoff = low_end
    if low_cutoff < 1:
        low_cutoff = 1
    if adata.shape[0] >= 40000 and low_cutoff <= 0.5:
        low_cutoff = 0.5
    num_low = int(np.ceil(np.sqrt(adata.shape[0]) * \
                          low_cutoff))
    eigvals_l, eigvecs_l = ss.linalg.eigsh(lap_mtx.astype(float),
                                           k=int(np.ceil(high_end * np.sqrt(adata.shape[0]))),
                                           which='LM')
    high_cutoff = np.ceil(kneed_select_values(eigvals_l, increasing=False) / \
                          np.sqrt(adata.shape[0]) * 1000) / 1000
    if high_cutoff < 1:
        high_cutoff = 1
    if high_cutoff >= high_end:
        high_cutoff = high_end
    if adata.shape[0] >= 40000 and high_cutoff <= 0.5:
        high_cutoff = 0.5
    num_high = int(np.ceil(np.sqrt(adata.shape[0]) * \
                           high_cutoff))

    adata.uns['FMs_after_select'] = {'low_FMs_frequency': eigvals_s[:num_low],
                                     'low_FMs': eigvecs_s[:, :num_low],
                                     'high_FMs_frequency': eigvals_l[(len(eigvals_l) - num_high):],
                                     'high_FMs': eigvecs_l[:, (len(eigvals_l) - num_high):]}

    return low_cutoff, high_cutoff


def detect_svg(adata,
               ratio_low_freq='infer',
               ratio_high_freq='infer',
               ratio_neighbors='infer',
               spatial_info=['array_row', 'array_col'],
               normalize_lap=False,
               filter_peaks=True,
               S=6,
               cal_pval=True):
    """
    Rank genes according to GFT score to find spatially variable genes based on
    graph Fourier transform.

    Parameters
    ----------
    adata : AnnData
        adata.X is the normalized count matrix. Besides, the spatial coordinates could be found in adata.obs or
        adata.obsm.
    ratio_low_freq : float | "infer", optional
        The ratio_low_freq will be used to determine the number of the FMs with
        low frequencies. Indeed, the ratio_low_freq * sqrt(number of spots) low
        frequency FMs will be calculated. If 'infer', the ratio_low_freq will be
        set to 1.0. The default is 'infer'.
    ratio_high_freq: float | 'infer', optional
        The ratio_high_freq will be used to determine the number of the FMs of
        high frequencies. Indeed, the ratio_high_freq * sqrt(number of spots) 
        high frequency FMs will be calculated. If 'infer', the ratio_high_freq
        will be set to 1.0. The default is 'infer'.
    ratio_neighbors: float | 'infer', optional
        The ratio_neighbors will be used to determine the number of neighbors
        when construct the KNN graph by spatial coordinates. Indeed, ratio_neighbors * sqrt(number of spots) / 2
         indicates the K. If 'infer', the parameter will be set to 1.0. The default is 'infer'.
    spatial_info : list | tuple | string, optional
        The column names of spatial coordinates in adata.obs_names or key
        in adata.varm_keys() to obtain spatial information. The default
        is ['array_row', 'array_col'].
    normalize_lap : bool, optional
        Whether you need to normalize laplacian matrix. The default is false.
    filter_peaks: bool, optional
        For calculated vectors/signals in frequency/spectral domain, whether
        filter low peaks to stress the important peaks. The default is True.
    S: int, optional
        The sensitivity parameter in Kneedle algorithm. A large S will enable
        more genes identified as svgs according to gft_score. The default is 6.
    cal_pval : bool, optional
        Whether you need to calculate p val by mannwhitneyu. The default is False.
    Returns
    -------
    score_df : dataframe
        Return gene information.

    """
    # Ensure parameters
    if ratio_low_freq == 'infer':
        num_low_frequency = int(np.ceil(np.sqrt(adata.shape[0])))
    else:
        num_low_frequency = int(np.ceil(np.sqrt(adata.shape[0]) * ratio_low_freq))
    if ratio_high_freq == 'infer':
        num_high_frequency = int(np.ceil(np.sqrt(adata.shape[0])))
    else:
        num_high_frequency = int(np.ceil(np.sqrt(adata.shape[0]) * ratio_high_freq))

    if ratio_neighbors == 'infer':
        num_neighbors = int(np.ceil(np.sqrt(adata.shape[0]) / 2))
    else:
        num_neighbors = int(np.ceil(np.sqrt(adata.shape[0]) / 2 * ratio_neighbors))
    if adata.shape[0] <= 500:
        num_neighbors = 4

    # Ensure gene index uniquely and all genes have expression  
    adata.var_names_make_unique()
    sc.pp.filter_genes(adata, min_cells=1)

    # Check dimensions
    if 'FMs_after_select' in adata.uns_keys():
        low_condition = (num_low_frequency == adata.uns['FMs_after_select']['low_FMs_frequency'].size)
        high_condition = (num_high_frequency == adata.uns['FMs_after_select']['high_FMs_frequency'].size)
    else:
        low_condition = False
        high_condition = False
    # ************ Construct graph and corresponding matrix *************
    lap_mtx = get_laplacian_mtx(adata, num_neighbors=num_neighbors,
                                spatial_key=spatial_info,
                                normalization=normalize_lap)

    # Next, calculate the eigenvalues and eigenvectors of the Laplacian
    # matrix as the Fourier modes with certain frequencies
    if not low_condition:
        eigvals_s, eigvecs_s = ss.linalg.eigsh(lap_mtx.astype(float),
                                               k=num_low_frequency,
                                               which='SM')
    else:
        eigvals_s, eigvecs_s = adata.uns['FMs_after_select']['low_FMs_frequency'], \
            adata.uns['FMs_after_select']['low_FMs']
        print('The precalculated low-frequency FMs are USED')
    if not high_condition:
        if num_high_frequency > 0:
            # Fourier bases of high frequency
            eigvals_l, eigvecs_l = ss.linalg.eigsh(lap_mtx.astype(float),
                                                   k=num_high_frequency,
                                                   which='LM')
    else:
        eigvals_l, eigvecs_l = adata.uns['FMs_after_select']['high_FMs_frequency'], \
            adata.uns['FMs_after_select']['high_FMs']
        print('The precalculated high-frequency FMs are USED')
    if num_high_frequency > 0:
        # eigenvalues
        eigvals = np.concatenate((eigvals_s, eigvals_l))
        # eigenvectors
        eigvecs = np.concatenate((eigvecs_s, eigvecs_l), axis=1)
    else:
        eigvals = eigvals_s
        eigvecs = eigvecs_s

    # ************************ Graph Fourier Transform *************************
    # Calculate GFT
    eigvecs_t = eigvecs.transpose()
    if type(adata.X) == np.ndarray:
        exp_mtx = preprocessing.scale(adata.X)
    else:
        exp_mtx = preprocessing.scale(adata.X.toarray())

    frequency_array = np.matmul(eigvecs_t, exp_mtx)
    frequency_array = np.abs(frequency_array)

    # Filter noise peaks
    if filter_peaks:
        frequency_array_thres_low = np.quantile(frequency_array[:num_low_frequency, :], q=0.5, axis=0)
        frequency_array_thres_high = np.quantile(frequency_array[num_low_frequency:, :], q=0.5, axis=0)
        for j in range(frequency_array.shape[1]):
            frequency_array[:num_low_frequency, :][frequency_array[:num_low_frequency, j] <= \
                                                   frequency_array_thres_low[j], j] = 0
            frequency_array[num_low_frequency:, :][frequency_array[num_low_frequency:, j] <= \
                                                   frequency_array_thres_high[j], j] = 0

    frequency_array = preprocessing.normalize(frequency_array,
                                              norm='l1',
                                              axis=0)

    eigvals = np.abs(eigvals)
    eigvals_weight = np.exp(-1 * eigvals)
    score_list = np.matmul(eigvals_weight, frequency_array)
    score_ave = np.matmul(eigvals_weight, (1 / len(eigvals)) * \
                          np.ones(len(eigvals)))
    score_list = score_list / score_ave
    print("Graph Fourier Transform finished!")

    # Rank genes according to smooth score
    adata.var["gft_score"] = score_list
    score_df = adata.var["gft_score"]
    score_df = pd.DataFrame(score_df)
    score_df = score_df.sort_values(by="gft_score", ascending=False)
    score_df.loc[:, "svg_rank"] = range(1, score_df.shape[0] + 1)
    adata.var["svg_rank"] = score_df.reindex(adata.var_names).loc[:, "svg_rank"]
    print("svg ranking could be found in adata.var['svg_rank']")

    # Determine cutoff of gft_score
    from kneed import KneeLocator
    magic = KneeLocator(score_df.svg_rank.values,
                        score_df.gft_score.values,
                        direction='decreasing',
                        curve='convex',
                        S=S)
    score_df['cutoff_gft_score'] = False
    score_df['cutoff_gft_score'][:(magic.elbow + 1)] = True
    adata.var['cutoff_gft_score'] = score_df['cutoff_gft_score']
    print("""The spatially variable genes judged by gft_score could be found 
          in adata.var['cutoff_gft_score']""")
    adata.varm['freq_domain_svg'] = frequency_array.transpose()
    print("""Gene signals in frequency domain when detect svgs could be found
          in adata.varm['freq_domain_svg']""")
    adata.uns["identify_svg_data"] = {}
    adata.uns["identify_svg_data"]['frequencies_low'] = eigvals_s
    adata.uns["identify_svg_data"]['frequencies_high'] = eigvals_l
    adata.uns["identify_svg_data"]['fms_low'] = eigvecs_s
    adata.uns["identify_svg_data"]['fms_high'] = eigvecs_l

    if cal_pval:
        if num_high_frequency == 0:
            raise ValueError("ratio_high_freq should be greater than 0")
        pval_list = test_significant_freq(
            freq_array=adata.varm['freq_domain_svg'],
            cutoff=num_low_frequency)
        from statsmodels.stats.multitest import multipletests
        qval_list = multipletests(np.array(pval_list), method='fdr_by')[1]
        adata.var['pvalue'] = pval_list
        adata.var['qvalue'] = qval_list
        score_df = adata.var.loc[score_df.index, :].copy()

    return score_df


def calculate_frequency_domain(adata,
                               ratio_low_freq='infer',
                               ratio_high_freq='infer',
                               ratio_neighbors='infer',
                               spatial_info=['array_row', 'array_col'],
                               return_freq_domain=True,
                               normalize_lap=False,
                               filter_peaks=False):
    """
    Obtain gene signals in frequency/spectral domain for all genes in 
    adata.var_names.

    Parameters
    ----------
    adata : AnnData
        adata.X is the normalized count matrix. Besides, the spatial coordinates could be found in adata.obs
        or adata.obsm.
    ratio_low_freq : float | "infer", optional
        The ratio_low_freq will be used to determine the number of the FMs with low frequencies. Indeed,
        the ratio_low_freq * sqrt(number of spots) low frequency FMs will be calculated. If 'infer',
        the ratio_low_freq will be set to 1.0. The default is 'infer'.
    ratio_high_freq: float | 'infer', optional
        The ratio_high_freq will be used to determine the number of the FMs with high frequencies. Indeed,
        the ratio_high_freq * sqrt(number of spots) high frequency FMs will be calculated. If 'infer',
        the ratio_high_freq will be set to 0. The default is 'infer'.
    ratio_neighbors: float | 'infer', optional
        The ratio_neighbors will be used to determine the number of neighbors when construct the KNN graph by spatial
        coordinates. Indeed, ratio_neighbors * sqrt(number of spots) / 2 indicates the K. If 'infer', the parameter
        will be set to 1.0. The default is 'infer'.
    spatial_info : list | tuple | str, optional
        The column names of spatial coordinates in adata.obs_keys() or a key in adata.obsm_keys.
        The default is ['array_row','array_col'].
    return_freq_domain : bool, optional
        Whether you need to return gene signals in frequency domain. The default is True.
    normalize_lap : bool, optional
        Whether you need to normalize laplacian matrix. The default is false.
    filter_peaks: bool, optional
        For calculated vectors/signals in frequency/spectral domain, whether filter low peaks to stress the important
        peaks. The default is False.

    Returns
    -------
    If return_freq_domain, return DataFrame, the index indicates the gene and 
    the columns indicates corresponding frequencies/smoothness.

    """
    # Critical parameters
    # Ensure parameters
    if ratio_low_freq == 'infer':
        num_low_frequency = int(np.ceil(np.sqrt(adata.shape[0])))
    else:
        num_low_frequency = int(np.ceil(np.sqrt(adata.shape[0]) * ratio_low_freq))
    if ratio_high_freq == 'infer':
        num_high_frequency = int(np.ceil(np.sqrt(adata.shape[0])))
    else:
        num_high_frequency = int(np.ceil(np.sqrt(adata.shape[0]) * ratio_high_freq))
    if adata.shape[0] >= 10000:
        num_high_frequency = int(np.ceil(np.sqrt(adata.shape[0])))

    if ratio_neighbors == 'infer':
        num_neighbors = int(np.ceil(np.sqrt(adata.shape[0]) / 2))
    else:
        num_neighbors = int(np.ceil(np.sqrt(adata.shape[0]) / 2 * ratio_neighbors))
    if adata.shape[0] <= 500:
        num_neighbors = 4

    # Ensure gene index uniquely and all gene had expression
    adata.var_names_make_unique()
    sc.pp.filter_genes(adata, min_cells=1)

    # ************** Construct graph and corresponding matrix ****************
    lap_mtx = get_laplacian_mtx(adata,
                                num_neighbors=num_neighbors,
                                spatial_key=spatial_info,
                                normalization=normalize_lap)

    # Calculate the eigenvalues and eigenvectors of the Laplace matrix
    np.random.seed(123)
    if num_high_frequency > 0:
        # Fourier modes of low frequency
        eigvals_s, eigvecs_s = ss.linalg.eigsh(lap_mtx.astype(float),
                                               k=num_low_frequency,
                                               which='SM')
        # Fourier mode of high frequency
        eigvals_l, eigvecs_l = ss.linalg.eigsh(lap_mtx.astype(float),
                                               k=num_high_frequency,
                                               which='LM')
        eigvals = np.concatenate((eigvals_s, eigvals_l))  # eigenvalues
        eigvecs = np.concatenate((eigvecs_s, eigvecs_l), axis=1)  # eigenvectors
    else:
        eigvals_s, eigvecs_s = ss.linalg.eigsh(lap_mtx.astype(float),
                                               k=num_low_frequency,
                                               which='SM')
        eigvecs = eigvecs_s
        eigvals = eigvals_s

    # ************************Graph Fourier Transform***************************
    # Calculate GFT
    eigvecs = eigvecs.transpose()
    if not ss.issparse(adata.X):
        exp_mtx = adata.X.copy()
    else:
        exp_mtx = adata.X.toarray().copy()
    exp_mtx = preprocessing.scale(exp_mtx, axis=0)
    frequency_array = np.matmul(eigvecs, exp_mtx)
    # Filter noise peaks
    if filter_peaks:
        frequency_array_thres = np.median(frequency_array, axis=0)
        for j in range(adata.shape[1]):
            frequency_array[frequency_array[:, j] <= \
                            frequency_array_thres[j], j] = 0
    # Spectral domain normalization
    frequency_array = preprocessing.normalize(frequency_array,
                                              norm='l1', axis=0)

    # ********************** Results of GFT ***********************************
    frequency_df = pd.DataFrame(frequency_array, columns=adata.var_names,
                                index=['low_spec_' + str(low) \
                                       for low in range(1, num_low_frequency + 1)] \
                                      + ['high_spec_' + str(high) \
                                         for high in range(1, num_high_frequency + 1)])
    adata.varm['freq_domain'] = frequency_df.transpose()
    adata.uns['frequencies'] = eigvals

    if return_freq_domain:
        return frequency_df


def freq2umap(adata,
              ratio_low_freq='infer',
              ratio_high_freq='infer',
              ratio_neighbors='infer',
              spatial_info=['array_row', 'array_col'],
              normalize_lap=False,
              filter_peaks=False):
    """
    Obtain gene signals in frequency/spectral domain for all genes in 
    adata.var_names and reduce dimension to 2 by UMAP.

    Parameters
    ----------
    adata : AnnData
        adata.X is the normalized count matrix. Besides, the spatial coordinates could be found in adata.obs
        or adata.obsm.
    ratio_low_freq : float | "infer", optional
        The ratio_low_freq will be used to determine the number of the FMs with low frequencies. Indeed,
        the ratio_low_freq * sqrt(number of spots) low frequency FMs will be calculated. If 'infer',
        the ratio_low_freq will be set to 1.0. The default is 'infer'.
    ratio_high_freq: float | 'infer', optional
        The ratio_high_freq will be used to determine the number of the FMs with high frequencies. Indeed,
        the ratio_high_freq * sqrt(number of spots) high frequency FMs will be calculated. If 'infer',
        the ratio_high_freq will be set to 0. The default is 'infer'.
    ratio_neighbors: float | 'infer', optional
        The ratio_neighbors will be used to determine the number of neighbors when construct the KNN graph by spatial
        coordinates. Indeed, ratio_neighbors * sqrt(number of spots) / 2 indicates the K. If 'infer', the parameter
        will be set to 1.0. The default is 'infer'.
    spatial_info : list | tuple | str, optional
        The column names of spatial coordinates in adata.obs_keys() or in adata.obsm_keys.
        The default is ['array_row','array_col'].
    normalize_lap : bool, optional
        Whether you need to normalize laplacian matrix. The default is false.
    filter_peaks: bool, optional
        For calculated vectors/signals in frequency/spectral domain, whether filter low peaks to stress the
        important peaks. The default is False.

    """
    if 'svg_rank' not in adata.var.columns:
        assert KeyError("adata.var['svg_rank'] is not available. Please run SpaGFT.rank_gene_smooth(adata) firstly.")
    ftup_adata = adata.copy()
    if 'log1p' in adata.uns_keys():
        ftup_adata.uns.pop('log1p')
    ftup_adata.X = adata.raw[:, adata.var_names].X
    sc.pp.log1p(ftup_adata)
    calculate_frequency_domain(ftup_adata,
                               ratio_low_freq=ratio_low_freq,
                               ratio_high_freq=ratio_high_freq,
                               ratio_neighbors=ratio_neighbors,
                               spatial_info=spatial_info,
                               return_freq_domain=False,
                               normalize_lap=normalize_lap,
                               filter_peaks=filter_peaks)
    adata.varm['gft_umap_svg'] = ftup_adata.varm['gft_umap']


def identify_ftu(adata,
                 svg_list='infer',
                 ratio_fms='infer',
                 ratio_neighbors=2,
                 spatial_info=['array_row', 'array_col'],
                 n_neighbors=15,
                 resolution=1,
                 weight_by_freq=False,
                 normalize_lap=False,
                 random_state=0,
                 **kwargs):
    """
    After identifying spatially variable genes, this function will group these spatially variable genes sharing common
    spatial patterns.

    Parameters
    ----------
    adata : AnnData
        adata.X is the normalized count matrix. Besides, the spatial coordinates could be found in adata.obs
        or adata.obsm.
    svg_list : list, optional
        The genes in svg_list will be grouped based on spatial patterns. The default is 'infer'.
    ratio_fms : float, optional
        The ratio_fms will be used to determine the number of the FMs with low frequencies. Indeed,
        the ratio_fms * sqrt(number of spots) low frequency FMs will be calculated. If 'infer',
        the ratio_low_freq will be determined automatically. The default is 'infer'.
    ratio_neighbors: float | 'infer', optional
        The ratio_neighbors will be used to determine the number of neighbors when construct the KNN graph by spatial
        coordinates. Indeed, ratio_neighbors * sqrt(number of spots) / 2 indicates the K. If 'infer', the parameter
        will be set to 2.0. The default is 2.
    spatial_info : list | tuple | str, optional
        The column names of spatial coordinates in adata.obs_keys() or a key in adata.obsm_keys.
         The default is ['array_row','array_col'].
    n_neighbors : int, optional
        The neighbors in gene similarity graph to perform louvain algorithm. 
        The default is 15.
    resolution : float | list | tuple, optional
        The resolution parameter in louvain algorithm. If resolution is float, resolution will be used directly.
        If resolution is a list, each value in this list will be used and the best value will be determined
        automatically. If resolution is tuple, it should be (start, end, step), and it is similar to a list.
        The default is 1.
    weight_by_freq : bool, optional
        Whether you need to weight FC according to frequencies. The default is False.
    normalize_lap : bool, optional
        Whether you need to normalize laplacian matrix. The default is false.
    random_state : int, optional
        The randomstate. The default is 0.
    **kwargs : kwargs 
        The parameters used in louvain algorithms and user can seek help in  sc.tl.louvain.

    Returns
    -------
    DataFrame
        The ftu information after identification.

    """
    # Find ftu by grouping Spatially variable genes with similar 
    # spatial patterns according to louvain algorithm.
    # Check conditions and determine parameters
    if isinstance(resolution, float) or isinstance(resolution, int):
        single_resolution = True
    else:
        single_resolution = False
    if isinstance(resolution, tuple) and len(resolution) == 3:
        start, end, step = resolution
        resolution = np.arange(start, end, step).tolist()
    elif isinstance(resolution, tuple):
        raise ValueError("""when resolution is a tuple, it should be (start, end, step)""")
    if isinstance(resolution, np.ndarray):
        resolution = resolution.tolist()

    assert isinstance(resolution, float) or isinstance(resolution, int) \
           or isinstance(resolution, list), 'please input resolution with type of float, int, list or tuple'

    if 'svg_rank' not in adata.var.columns:
        assert KeyError("adata.var['svg_rank'] is not available. Please run SpaGFT.detect_svg(adata) firstly.")
    if ratio_fms == 'infer':
        if adata.shape[0] <= 500:
            ratio_fms = 4
        elif adata.shape[0] <= 10000:
            ratio_fms = 2
        else:
            ratio_fms = 1
    if isinstance(svg_list, str):
        if svg_list == 'infer':
            gene_score = adata.var.sort_values(by='svg_rank')
            adata = adata[:, gene_score.index]
            svg_list = adata.var[adata.var.cutoff_gft_score][adata.var.qvalue < 0.05].index.tolist()
    ftup_adata = adata[:, svg_list].copy()
    if 'log1p' in adata.uns_keys():
        ftup_adata.uns.pop('log1p')
    ftup_adata.X = adata[:, svg_list].X.copy()
    calculate_frequency_domain(ftup_adata,
                               ratio_low_freq=ratio_fms,
                               ratio_high_freq=0,
                               ratio_neighbors=ratio_neighbors,
                               spatial_info=spatial_info,
                               return_freq_domain=False,
                               normalize_lap=normalize_lap,
                               filter_peaks=False)
    # Create new anndata to store freq domain information
    gft_adata = sc.AnnData(ftup_adata.varm['freq_domain'])
    if weight_by_freq:
        weight_list = 1 / (1 + 0.01 * ftup_adata.uns['frequencies'])
        gft_adata.X = np.multiply(gft_adata.X, weight_list)
        gft_adata.X = preprocessing.normalize(gft_adata.X, norm='l1')
    # clustering
    gft_adata = gft_adata[svg_list, :]
    adata.uns['detect_ftu_data'] = {}
    adata.uns['detect_ftu_data']['freq_domain_svgs'] = ftup_adata[:, svg_list].varm['freq_domain']
    if gft_adata.shape[1] >= 400:
        sc.pp.pca(adata)
        sc.pp.neighbors(gft_adata, n_neighbors=n_neighbors)
    else:
        sc.pp.neighbors(gft_adata, n_neighbors=n_neighbors, use_rep='X')

    # Determining the resolution data type, if resolution is the type of list, 
    # we will select the optimal resolution
    if isinstance(resolution, list):
        # The minimum value of cosine similarity of overlap: 
        # count_ftus_select_cs_score
        count_ftus_select_cs_score = 1
        best_resolution = None
        overlap_scores = {}
        # Iterate through the list of resolution and calculate the 
        # clustering results
        for resolution_index, resolution_value in enumerate(resolution):
            gft_adata_current = gft_adata.copy()
            sc.tl.louvain(gft_adata_current,
                          resolution=resolution_value,
                          random_state=random_state,
                          key_added='louvain',
                          **kwargs)

            gft_adata_current.obs.louvain = [str(eval(i_ftu) + 1) for i_ftu in gft_adata_current.obs.louvain.tolist()]
            gft_adata_current.obs.louvain = pd.Categorical(gft_adata_current.obs.louvain)

            # ftu pseudo expression
            all_ftus_current = gft_adata_current.obs.louvain.cat.categories
            ftu_df_current = pd.DataFrame(0,
                                          index=ftup_adata.obs_names,
                                          columns='ftu_' + all_ftus_current)
            # Calculate the clustering of each ftu
            for ftu in all_ftus_current:
                pseudo_exp = ftup_adata[:,
                             gft_adata_current.obs.louvain[gft_adata_current.obs.louvain == ftu].index].X.sum(axis=1)
                pseudo_exp = np.ravel(pseudo_exp)
                # Calculate the clustering results
                predict_ftu = KMeans(n_clusters=2, random_state=random_state).fit_predict(pseudo_exp.reshape(-1, 1))

                # Correct clustering results
                pseudo_exp_median = np.median(pseudo_exp)
                pseudo_exp_cluster = np.where(pseudo_exp > pseudo_exp_median, 1, 0)

                cluster_middle_param = sum(abs(predict_ftu - pseudo_exp_cluster))
                cluster_middle_param_reverse = sum(abs(predict_ftu - abs(pseudo_exp_cluster - 1)))
                if cluster_middle_param > cluster_middle_param_reverse:
                    predict_ftu = abs(predict_ftu - 1)
                ftu_df_current['ftu_' + str(ftu)] = predict_ftu

            # Correct cosine similarity of overlap for clustering results
            overlap_cs_score = get_overlap_cs_core(ftu_df_current.values.T)
            print("""resolution: %.3f;  """ % resolution_value + """score: %.4f""" % overlap_cs_score)
            overlap_scores['res_' + '%.3f' % resolution_value] = np.round(overlap_cs_score * 1e5) / 1e5
            # select the optimal resolution
            if count_ftus_select_cs_score > overlap_cs_score or resolution_index == 0:
                count_ftus_select_cs_score = overlap_cs_score
                best_resolution = resolution_value
            resolution = best_resolution

    # Next, clustering genes for given resolution
    sc.tl.louvain(gft_adata, resolution=resolution, random_state=random_state,
                  **kwargs, key_added='louvain')
    sc.tl.umap(gft_adata)
    adata.uns['detect_ftu_data']['gft_umap_ftu'] = pd.DataFrame(gft_adata.obsm['X_umap'],
                                                                index=gft_adata.obs.index,
                                                                columns=['UMAP_1', 'UMAP_2'])
    gft_adata.uns['gft_genes_ftu'] = [str(eval(i_ftu) + 1) for i_ftu in gft_adata.obs.louvain.tolist()]
    gft_adata.obs.louvain = [str(eval(i_ftu) + 1) for i_ftu in gft_adata.obs.louvain.tolist()]
    gft_adata.obs.louvain = pd.Categorical(gft_adata.obs.louvain)
    adata.var['ftu'] = 'None'
    adata.var.loc[gft_adata.obs_names, 'ftu'] = gft_adata.obs.louvain
    adata.var['ftu'] = pd.Categorical(adata.var['ftu'])

    # ftu pseudo expression
    all_ftus = gft_adata.obs.louvain.cat.categories
    ftu_df = pd.DataFrame(0, index=adata.obs_names, columns='ftu_' + all_ftus)
    pseudo_df = pd.DataFrame(0, index=adata.obs_names, columns='ftu_' + all_ftus)
    for ftu in all_ftus:
        pseudo_exp = ftup_adata[:,
                     gft_adata.obs.louvain[gft_adata.obs.louvain == ftu].index].X.sum(axis=1)
        pseudo_exp = np.ravel(pseudo_exp)
        pseudo_df['ftu_' + str(ftu)] = pseudo_exp.copy()
        predict_ftu = KMeans(n_clusters=2, random_state=random_state).fit_predict(pseudo_exp.reshape(-1, 1))

        pseudo_exp_median = np.median(pseudo_exp)
        pseudo_exp_cluster = np.where(pseudo_exp > pseudo_exp_median, 1, 0)

        cluster_middle_param = sum(abs(predict_ftu - pseudo_exp_cluster))
        cluster_middle_param_reverse = sum(abs(predict_ftu - abs(pseudo_exp_cluster - 1)))
        if cluster_middle_param > cluster_middle_param_reverse:
            predict_ftu = abs(predict_ftu - 1)
        ftu_df['ftu_' + str(ftu)] = predict_ftu
        adata.obsm['ftu_pseudo_expression'] = pseudo_df.copy()

    ftu_df = ftu_df.astype(str)
    adata.obsm['ftu_binary'] = ftu_df.copy()
    # Obtain freq signal
    freq_signal_ftu_df = pd.DataFrame(0, index=ftu_df.columns,
                                      columns=ftup_adata.varm['freq_domain'].columns)

    for ftu in all_ftus:
        ftu_gene_list = gft_adata.obs.louvain[gft_adata.obs.louvain == ftu].index
        freq_signal = ftup_adata.varm['freq_domain'].loc[ftu_gene_list,
                      :].sum(axis=0)
        freq_signal = np.abs(freq_signal)
        freq_signal = freq_signal / sum(freq_signal)
        freq_signal_ftu_df.loc['ftu_' + ftu, :] = freq_signal
    adata.uns['detect_ftu_data']['freq_signal_ftu'] = freq_signal_ftu_df
    adata.uns['detect_ftu_data']['low_freq_domain_svg'] = pd.DataFrame(gft_adata.X.copy(),
                                                                       index=gft_adata.obs_names,
                                                                       columns=['low_freq_' + str(i + 1) \
                                                                                for i in range(gft_adata.shape[1])])
    if not single_resolution:
        adata.uns['detect_ftu_data']['overlap_curve'] = pd.DataFrame(overlap_scores, index=['score'])
        adata.uns['detect_ftu_data']['overlap_curve'] = adata.uns['detect_ftu_data']['overlap_curve'].transpose()

    return adata.var.loc[svg_list, :].copy(), adata
