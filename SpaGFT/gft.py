#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb  9 18:47:29 2022

@author: frank
"""
import pandas as pd
import numpy as np
import warnings
import scipy.sparse as ss
import scanpy as sc
from sklearn import preprocessing
from SpaGFT.utilis import get_laplacian_mtx
warnings.filterwarnings("ignore")


def permutation_signal(signal_array, num_permutaion=1000):
    """
    Permutate gene score in spatial domain randomly.

    Parameters
    ----------
    signal_array : list | array
        A one-dimensional array indicate gene expression on all spots.
    num_permutaion : int, optional
        The number of permutation. The default is 1000.

    Returns
    -------
    total_signals : array
        The permutaed gene expression signals.

    """
    signal_array = np.array(signal_array)
    total_signals = signal_array * np.ones((num_permutaion,
                                             len(signal_array)))
    
    for i in range(num_permutaion):
        total_signals[i, :] = np.random.permutation(total_signals[i, :])
        
    return total_signals

# def signficant_test(exp_mtx, gene_score, eigvals, eigvecs_T, num_permutaion=1000,
#                     num_pool=100, spec_norm='l1'):
#     from multiprocessing import Pool
#     from multiprocessing.dummy import Pool as ThreadPool
#     from scipy.stats import ttest_1samp

#     def test_by_permutaion(gene_index):
#         graph_signal = exp_mtx[gene_index, :]
#         total_signals = permutation_signal(signal_array=graph_signal,
#                                            num_permutaion=num_permutaion)
#         frequency_array = np.matmul(eigvecs_T, total_signals.transpose())
#         frequency_array = np.abs(frequency_array)
#         if spec_norm != None:
#             frequency_array = preprocessing.normalize(frequency_array, 
#                                                       norm=spec_norm,
#                                                       axis=0)
#         score_list = np.matmul(2 ** (-1 * eigvals), frequency_array)
#         # score_max = 2 ** (-1 * eigvals[0])       
#         score_list = score_list / score_max
#         pval = ttest_1samp(score_list, gene_score[gene_index], alternative='less').pvalue
#         return pval
        
#     score_max = np.matmul(2 ** (-1 * eigvals), (1/len(eigvals)) * \
#                           np.ones(len(eigvals)))  
#     gene_index_list = list(range(exp_mtx.shape[0]))
#     pool = ThreadPool(num_pool)
#     res = pool.map(test_by_permutaion, gene_index_list)
    
#     return res

def my_eigsh(args_tupple):
    """
    The function is used to multi-process calculate using Pool.

    Parameters
    ----------
    args_tupple : tupple 
        The args_tupple contains three elements, that are, Lplacian matrix, k 
        and which.

    Returns
    -------
    (eigvals, eigvecs)

    """
    lap_mtx=args_tupple[0]
    k=args_tupple[1]
    which=args_tupple[2]
    eigvals, eigvecs = ss.linalg.eigsh(lap_mtx.astype(float), 
                                           k=k,
                                           which=which) 
    return((eigvals, eigvecs))

def low_frequency_filter(adata,
                         num_low_frequency='infer',
                         num_high_frequency='infer',
                         c = 0.001,
                         num_neighbors='infer',
                         spatial_names = ['array_row', 'array_col'],
                         spatial_key=None,
                         normalization=False,
                         inplace=False):
    """
    Implement gene expression with low frequency filter. After this step, the 
    spatially variables genes will be more smooth than previous. The functions
    can also be treated as denoising. Note that the denosing results is related
    to spatial graph topology so that only the resulsts of spatially variable
    genes could be convincing.

    Parameters
    ----------
    adata : AnnData
        adata.X is the normalized count matrix. Besides, the spatial coordinat-
        es could be found in adata.obs or adata.obsm.
    num_low_frequency : int, optional
        The number of low frequencys. The default is 100.
    num_high_frequency : int, optional
        The number of high frequencys. The default is 100.
    c: float, optional
        c balance the smoothness and difference with previous expresssion.
        A high can achieve better smothness. c should be setted to [0, 0.05].
    num_neighbors: int, optional
        The number of neighbors for each node/spot/pixel when contrcut graph. 
        The defalut if 6.
    spatial_name : list or tupple, optional
        the column names of spaital coordinates in adata.obs_names. The default
        is ['array_row', 'array_col'].
    spatial_key=None : None | string
        Get the coordinate information by adata.obsm[spaital_key]. If none, 
        the spatial information could be found by adata.obs[spatial_names].
    normlization : bool. optional
        Whether need to normalize the Laplcian matrix.
        

    Returns
    -------
    count_matrix: DataFrame

    """
    if num_low_frequency == 'infer':
        if adata.shape[0] <= 400:
            num_low_frequency = min(adata.shape[0], 300)
        elif adata.shape[0] < 1000:
            num_low_frequency = min(adata.shape[0], 800)
        elif adata.shape[0] < 5000:
            num_low_frequency = min(adata.shape[0], 1400)
        else:
            num_low_frequency = 2000
    if num_neighbors == 'infer':
        if adata.shape[0] <= 400:
            num_neighbors = 4
        elif adata.shape[0] < 1000:
            num_neighbors = 16
        elif adata.shape[0] < 5000:
            num_neighbors = 32
        else:
            num_neighbors = 128
    if num_high_frequency == 'infer':
        num_high_frequency = 0
        
    adata.var_names_make_unique()
    sc.pp.filter_genes(adata, min_cells=1)
    # Get Laplacian matrix according to coordinates 
    lap_mtx = get_laplacian_mtx(adata, num_neighbors=num_neighbors,
                                spatial_names=spatial_names,
                                spatial_key=spatial_key,
                                normalization=normalization)
    
    # Fourier bases of low frequency
    eigvals_s, eigvecs_s = ss.linalg.eigsh(lap_mtx.astype(float), 
                                           k=num_low_frequency,
                                           which='SM') 
    if num_high_frequency > 0:
        # Fourier bases of high frequency    
        eigvals_l, eigvecs_l = ss.linalg.eigsh(lap_mtx.astype(float),
                                               k=num_high_frequency, 
                                               which='LM')
        eigvals = np.concatenate((eigvals_s, eigvals_l))         # eigenvalues
        eigvecs = np.concatenate((eigvecs_s, eigvecs_l), axis=1) # eigenvectors
    else:
        eigvals = eigvals_s
        eigvecs = eigvecs_s
    
    # ************************Graph Fourier Tranform***************************
    # Calculate GFT
    eigvecs_T = eigvecs.transpose()
    if isinstance(adata.X, np.ndarray):
        exp_mtx = adata.X
    else:
        exp_mtx = adata.X.toarray()
    frequency_array = np.matmul(eigvecs_T, exp_mtx)
    filter_list = [1 / (1 + c * eigv) for eigv in eigvals]
    filter_array = np.matmul(np.diag(filter_list), frequency_array)
    filter_array = np.matmul(eigvecs, filter_array)
    if inplace:
        adata.X = filter_array
    filter_array = pd.DataFrame(filter_array, index=adata.obs_names,
                                columns=adata.var_names)
    return filter_array

def rank_gene_smooth(adata, num_low_frequency='infer', 
                     num_high_frequency='infer', 
                     num_neighbors='infer',
                     spatial_names=['array_row', 'array_col'],
                     spatial_key=None,
                     normalization_lap=False,
                     exp_norm='l2',
                     filter_peaks=False,
                     q=0.2,
                     spec_norm='l1', 
                     multi_process=False):
    """
    Rank genes to find spatially variable genes by graph Fourier transform.

    Parameters
    ----------
    adata : AnnData
        adata.X is the normalized count matrix. Besides, the spatial coordinat-
        es could be found in adata.obs or adata.obsm.
    num_low_frequency : int | 'infer', optional
        The number of low frequencys. The default is 'infer'.
    num_high_frequency : int | 'infer', optional
        The number of high frequencys. The default is 'infer'.
    num_neighbors: int | 'infer', optional
        The number of neighbors for each node/spot/pixel when contrcut graph. 
        The defalut if 'infer'.
    spatial_name : list or tupple, optional
        The column names of spaital coordinates in adata.obs_names. The default
        is ['array_row', 'array_col'].
    spatial_key : string or None,
        If spatial_keq is provided, the spatial information could be found in
        adata.obsm['spaital_key'] (spot * coordinate).
    normalization_lap : bool, optional
        Whether need to normalize laplacian matrix. The default is false.
    exp_norm : None or {'l1', 'l2', 'max'}, optional
        If None, the expression matrix will not be normzalized. Otherwise, the 
        given exp_norm method will be implemented.
    filter_peaks: bool, optional
        For calculated vectors/signals in frequency/spectral domian, whether
        filter low peaks to stress the important peaks. The default is False.
    q: float, optional
        Quantile or sequence of quantiles to compute in frequency/spectral 
        domian. The default is 0.5.
    spec_norm: None or {'l1', 'l2', 'max'}, optional
        If None, the signal in spectral domian matrix will not be normzalized. 
        Otherwise, the given spec_norm method will be implemented. The default
        is l1. 
    multi_process: bool
        Whether need to multiple processes when calculate the eigenvalues and 
        eigenvectors. The default is False.

    Returns
    -------
    score_df : dataframe
        Return gene information.

    """
    # Ensure parameters
    if num_low_frequency == 'infer':
        if adata.shape[0] <= 400:
            num_low_frequency = min(adata.shape[0], 50)
        elif adata.shape[0] < 1000:
            num_low_frequency = min(adata.shape[0], 70)
        elif adata.shape[0] < 5000:
            num_low_frequency = min(adata.shape[0], 100)
        elif adata.shape[0] < 100000:
            num_low_frequency = 400
        else:
            num_low_frequency = 400
    if num_neighbors == 'infer':
        if adata.shape[0] <= 400:
            num_neighbors = 4
        elif adata.shape[0] < 1000:
            num_neighbors = 16
        elif adata.shape[0] < 5000:
            num_neighbors = 32
        elif adata.shape[0] < 50000:
            num_neighbors = 128
        else:
            num_neighbors = 256
    if num_high_frequency == 'infer':
        num_high_frequency = 0
        
    
    adata.var_names_make_unique()
    sc.pp.filter_genes(adata, min_cells=1)
    # ***************Construct graph and corresponding matrixs*****************
    lap_mtx = get_laplacian_mtx(adata, num_neighbors=num_neighbors,
                                spatial_names=spatial_names,
                                spatial_key=spatial_key,
                                normalization=normalization_lap)
    
    
    # Next, calculate the eigenvalues and eigenvectors of the Laplace matrix
    np.random.seed(123) 
    if num_high_frequency > 0:
        if multi_process == True:
            from multiprocessing import Pool
            import gc
            pool = Pool(1)
            out = pool.map_async(my_eigsh, [(lap_mtx.astype(float), 
                                         num_high_frequency,
                                         'LM'),])
            pool.close()
            eigvals_s, eigvecs_s = ss.linalg.eigsh(lap_mtx.astype(float), 
                                                   k=num_low_frequency,
                                                   which='SM')  
            pool.join()
            eigvals_l, eigvecs_l = out.get()[0]
            del out
            gc.collect()
        else:
            # Fourier bases of low frequency
            eigvals_s, eigvecs_s = ss.linalg.eigsh(lap_mtx.astype(float), 
                                                   k=num_low_frequency,
                                                   which='SM')  
            # Fourier bases of high frequency
            eigvals_l, eigvecs_l = ss.linalg.eigsh(lap_mtx.astype(float),
                                                   k=num_high_frequency, 
                                                   which='LM')       
        eigvals = np.concatenate((eigvals_s, eigvals_l))         # eigenvalues
        eigvecs = np.concatenate((eigvecs_s, eigvecs_l), axis=1) # eigenvectors
    else:
        # Fourier bases of low frequency
        eigvals_s, eigvecs_s = ss.linalg.eigsh(lap_mtx.astype(float), 
                                               k=num_low_frequency,
                                               which='SM') 
        eigvals = eigvals_s
        eigvecs = eigvecs_s
    
    # ************************Graph Fourier Tranform***************************
    # Calculate GFT
    eigvecs_T = eigvecs.transpose()
    if exp_norm != None:
        # exp_mtx = preprocessing.normalize(adata.X.toarray(), axis=0, 
        #                                   norm=exp_norm)
        if type(adata.X) == np.ndarray:
            exp_mtx = preprocessing.scale(adata.X)
        else:
            exp_mtx = preprocessing.scale(adata.X.toarray())
    frequency_array = np.matmul(eigvecs_T, exp_mtx)
    frequency_array = np.abs(frequency_array)
    
    # Filter noise peaks
    if filter_peaks == True:
        frequency_array_thres = np.quantile(frequency_array, q=q, axis=0)
        for j in range(adata.shape[1]):
            frequency_array[frequency_array[:, j] <= \
                            frequency_array_thres[j], j] = 0
    # spectral domian normalization
    if spec_norm != None:
        frequency_array = preprocessing.normalize(frequency_array, 
                                                  norm=spec_norm,
                                                  axis=0)
    score_list = np.matmul(2 ** (-1 * eigvals), frequency_array)
    # score_max = 2 ** (-1 * eigvals[0])  
    score_max = np.matmul(2 ** (-1 * eigvals), (1/len(eigvals)) * \
                          np.ones(len(eigvals)))       
    score_list = score_list / score_max
    print("Graph Fourier Transform finished!")
    
    # ****************** calculate pval ***************************
    # pval_list = signficant_test(exp_mtx=exp_mtx.transpose()[:400, :],
    #                             gene_score=score_list[:400], eigvals=eigvals,
    #                             eigvecs_T=eigvecs_T)
    # print(pval_list)

    # Rank genes according to smooth score
    adata.var["smooth_score"] = score_list 
    score_df =adata.var["smooth_score"]
    score_df =pd.DataFrame(score_df)
    score_df = score_df.sort_values(by="smooth_score", ascending=False) 
    
    score_df.loc[:, "SVG_Rank"] = range(1, score_df.shape[0] + 1)
    adata.var["SVG_Rank"] = score_df.reindex(adata.var_names).loc[:, "SVG_Rank"]
    print("SVG ranking could be found in adata.obs['Rank']")
    
    adata.varm['freq_domain'] = frequency_array.transpose()
    print("Gene signals in frequency domain could be found in \
          adata.varm['freq_domain']")
    
    return score_df


def calculate_frequcncy_domain(adata, num_low_frequency='infer', 
                               num_high_frequency='infer', 
                               num_neighbors='infer',
                               spatial_names=['array_row', 'array_col'],
                               spatial_key=None,
                               normalization_lap=False,
                               exp_norm='l2',
                               filter_peaks=True,
                               q=0.5,
                               spec_norm='l1',
                               multi_process=False):
    """
    Obtain gene signals in frequency/spectral domian fo given genes in 
    adata.var_names for given graph topology.

    Parameters
    ----------
    adata : AnnData
        adata.X is the normalized count matrix. Besides, the spatial coordinat-
        es could be found in adata.obs or adata.obsm.
    num_low_frequency : int | 'infer', optional
        The number of low frequencys. The default is 'infer'.
    num_high_frequency : int | 'infer', optional
        The number of high frequencys. The default is 'infer'.
    num_neighbors: int | 'infer', optional
        The number of neighbors for each node/spot/pixel when contrcut graph. 
        The defalut if 'infer'.
    spatial_name : list or tupple, optional
        The column names of spaital coordinates in adata.obs_names. The default
        is ['array_row', 'array_col'].
    spatial_key : string or None,
        If spatial_keq is provided, the spatial information could be found in
        adata.obsm['spaital_key'] (spot * coordinate).
    normalization_lap : bool, optional
        Whether need to normalize laplacian matrix. The default is false.
    exp_norm : None or {'l1', 'l2', 'max'}, optional
        If None, the expression matrix will not be normzalized. Otherwise, the 
        given exp_norm method will be implemented.
    filter_peaks: bool, optional
        For calculated vectors/signals in frequency/spectral domian, whether
        filter low peaks to stress the important peaks. The default is False.
    q: float, optional
        Quantile or sequence of quantiles to compute in frequency/spectral 
        domian. The default is 0.5.
    spec_norm: None or {'l1', 'l2', 'max'}, optional
        If None, the signal in spectral domian matrix will not be normzalized. 
        Otherwise, the given spec_norm method will be implemented. The default
        is l1. 
    multi_process: bool
        Whether need to multiple processes when calculate the eigenvalues and 
        eigenvectors. The default is False.

    Returns
    -------
    DataFrame, the index indicates the gene and the columns indicates correspo-
    nding frequecies/smoothness. 

    """ 
    # critical parameters
    if num_low_frequency == 'infer':
        if adata.shape[0] <= 400:
            num_low_frequency = min(adata.shape[0], 100)
        elif adata.shape[0] < 1000:
            num_low_frequency = 200
        elif adata.shape[0] < 5000:
            num_low_frequency = 400
        else:
            num_low_frequency = 1000
    if num_neighbors == 'infer':
        if adata.shape[0] <= 400:
            num_neighbors = 4
        elif adata.shape[0] < 1000:
            num_neighbors = 16
        elif adata.shape[0] < 5000:
            num_neighbors = 32
        else:
            num_neighbors = 128
    if num_high_frequency == 'infer':
        num_high_frequency = 0
       
    adata.var_names_make_unique()
    sc.pp.filter_genes(adata, min_cells=1)
    # ***************Construct graph and corresponding matrixs*****************
    lap_mtx = get_laplacian_mtx(adata, num_neighbors=num_neighbors,
                                spatial_names=spatial_names,
                                spatial_key=spatial_key,
                                normalization=normalization_lap)
    
    
    # Calculate the eigenvalues and eigenvectors of the Laplace matrix
    np.random.seed(123)
    if num_high_frequency > 0: 
        if multi_process == True:   # whether need to multi-process
            from multiprocessing import Pool
            import gc
            pool = Pool(1)
            out = pool.map_async(my_eigsh, [(lap_mtx.astype(float), 
                                         num_high_frequency,
                                         'LM'), ])
            pool.close()
            eigvals_s, eigvecs_s = ss.linalg.eigsh(lap_mtx.astype(float), 
                                                   k=num_low_frequency,
                                                   which='SM')  
            pool.join()
            eigvals_l, eigvecs_l = out.get()[0]
            del out
            gc.collect()
        else:
            # Fourier bases of low frequency
            eigvals_s, eigvecs_s = ss.linalg.eigsh(lap_mtx.astype(float), 
                                                   k=num_low_frequency,
                                                   which='SM')  
            # Fourier bases of high frequency    
            eigvals_l, eigvecs_l = ss.linalg.eigsh(lap_mtx.astype(float),
                                                   k=num_high_frequency, 
                                                   which='LM')
        # eigvals = np.concatenate((eigvals_s, eigvals_l))            # eigenvalues
        eigvecs = np.concatenate((eigvecs_s, eigvecs_l), axis=1)    # eigenvectors
    else:
        eigvals_s, eigvecs_s = ss.linalg.eigsh(lap_mtx.astype(float), 
                                               k=num_low_frequency,
                                               which='SM')
        eigvecs = eigvecs_s
   
    
    # ************************Graph Fourier Tranform***************************
    # Calculate GFT
    eigvecs = eigvecs.transpose()
    if exp_norm != None:
        if type(adata.X) == np.ndarray:
            exp_mtx = preprocessing.scale(adata.X)
        else:
            exp_mtx = preprocessing.scale(adata.X.toarray())
    frequency_array = np.matmul(eigvecs, exp_mtx)
    frequency_array = np.abs(frequency_array)
    
    # Filter noise peaks
    if filter_peaks == True:
        frequency_array_thres = np.quantile(frequency_array, q=q, axis=0)
        for j in range(adata.shape[0]):
            frequency_array[frequency_array[:, j] <= frequency_array_thres[j], j] = 0
    # spectral domian normalization
    if spec_norm != None:
        frequency_array = preprocessing.normalize(frequency_array, 
                                                  norm=spec_norm,
                                                  axis=0)
    
    # **********************results of GFT*************************************
    frequency_df = pd.DataFrame(frequency_array, columns=adata.var_names, 
                                index=['low_spec_' + str(low) \
                                         for low in range(1, num_low_frequency + 1)] \
                                    + ['high_spec_' + str(high) \
                                       for high in range(1, num_high_frequency + 1)])
    adata.varm['freq_domain'] = frequency_df.values.transpose()
    
    return frequency_df


    
    
    
    
    
    
    
    