import pandas as pd
import numpy as np
import warnings
import scipy.sparse as ss
import scanpy as sc
from sklearn import preprocessing
from SpaGFT.utilis import get_laplacian_mtx
warnings.filterwarnings("ignore")


def _correct_pvalues_for_multiple_testing(pvalues, 
                                          correction_type = "Benjamini-Hochberg"): 
    """
    Correct pvalues to obtain the adjusted pvalues 

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
        values = [ (pvalue, i) for i, pvalue in enumerate(pvalues) ]                                      
        values.sort()
        for rank, vals in enumerate(values):                                                              
            pvalue, i = vals
            new_pvalues[i] = (n-rank) * pvalue                                                            
    elif correction_type == "Benjamini-Hochberg":                                                         
        values = [ (pvalue, i) for i, pvalue in enumerate(pvalues) ]                                      
        values.sort()
        values.reverse()                                                                                  
        new_values = []
        for i, vals in enumerate(values):                                                                 
            rank = n - i
            pvalue, index = vals                                                                          
            new_values.append((n/rank) * pvalue)                                                          
        for i in range(0, int(n)-1):  
            if new_values[i] < new_values[i+1]:                                                           
                new_values[i+1] = new_values[i]                                                           
        for i, vals in enumerate(values):
            pvalue, index = vals
            new_pvalues[index] = new_values[i]                                                                                                                  
    return new_pvalues

def _permutation_signal(signal_array, num_permutaion=1000):
    """
    Permutate gene signals in spatial domain randomly.

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

def _significant_test_permutation(exp_mtx, 
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
        The count matrix of gene expresssion. (spots * genes)
    gene_score : 1D-array
        The calculated gene scores. 
    eigvals : array
        The eigenvalues of Laplacian matrix.
    eigvecs_T : array
        The eigenvectors of Laplacian matrix.
    num_permutaion : int, optional
        The number of permutations. The default is 1000.
    num_pool : int, optional
        The cores used for umltiprocess calculation to accelerate speed. The 
        default is 200.
    spec_norm : str, optional
        The method to normalize graph signals in spectral domain. The default 
        is 'l1'.

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
        total_signals = _permutation_signal(signal_array=graph_signal,
                                            num_permutaion=num_permutaion)
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
        
    score_max = np.matmul(2 ** (-2 * eigvals), (1/len(eigvals)) * \
                          np.ones(len(eigvals)))  
    gene_index_list = list(range(exp_mtx.shape[0]))
    pool = ThreadPool(num_pool)
    res = pool.map(_test_by_permutaion, gene_index_list)
    
    return res

def _test_significant_freq(freq_array,
                           cutoff,
                           num_pool=200,
                           combine_method='fisher'):
    """
    Significance test by camparing the intensities in low frequency FMs and 
    in high frequency FMs. 

    Parameters
    ----------
    freq_array : array
        The graph signals of genes in frequency domain. 
    cutoff : int
        Watershed between low frequency signals and high frequency signals.
    q : float, optional
        The quantile used in filtering noise peaks. The default is 0.3.
    num_pool : int, optional
        The cores used for umltiprocess calculation to accelerate speed. The 
        default is 200.
    combine_method : str, optional
        The method for intergrating pvalues from multi-tests. The default is 
        'fisher'.

    Returns
    -------
    array
        The calculated p values.

    """
    from scipy.stats import ranksums
    from multiprocessing.dummy import Pool as ThreadPool
    
    def _test_by_feq(gene_index):
        freq_signal = freq_array[gene_index, :]
        freq_1 = freq_signal[:cutoff]
        freq_1 = freq_1[freq_1 > 0]
        freq_2 = freq_signal[cutoff:]
        freq_2 = freq_2[freq_2 > 0]
        if freq_1.size <= 15 or freq_2.size <= 15:
            freq_1 = np.concatenate((freq_1, freq_1, freq_1, freq_1))
            freq_2 = np.concatenate((freq_2, freq_2, freq_2, freq_2))
        elif freq_1.size <= 50 or freq_2.size <= 50:
            freq_1 = np.concatenate((freq_1, freq_1))
            freq_2 = np.concatenate((freq_2, freq_2))
        pval = ranksums(freq_1, freq_2, alternative='greater').pvalue
        return pval
    
    gene_index_list = list(range(freq_array.shape[0]))
    pool = ThreadPool(num_pool)
    res = pool.map(_test_by_feq, gene_index_list)
    
    return res
    
def _my_eigsh(args_tupple):
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

def low_pass_imputation(adata,
                        ratio_low_freq='infer',
                        ratio_high_freq='infer',
                        ratio_neighbors='infer',
                        c=0.0001,
                        spatial_info = ['array_row', 'array_col'],
                        normalize_lap=False,
                        inplace=False):
    """
    Implement gene expression with low-pass filter. After this step, the 
    spatially variables genes will be more smooth than previous. The function
    can also be treated as denoising. Note that the denosing results is related
    to spatial graph topology so that only the resulsts of spatially variable
    genes could be convincing.

    Parameters
    ----------
    adata : AnnData
        adata.X is the normalized count matrix. Besides, the spatial coordinat-
        es of all spots could be found in adata.obs or adata.obsm.
    ratio_low_freq : float | "infer", optional
        The ratio_low_freq will be used to determine the number of the FMs with
        low frequencies. Indeed, the ratio_low_freq * sqrt(number of spots) low
        frequecy FMs will be calculated. If 'infer', the ratio_low_freq will be
        set to 1.0. The default is 'infer'.
    ratio_high_freq: float | 'infer', optional
        The ratio_high_freq will be used to determine the number of the FMs with
        high frequencies. Indeed, the ratio_high_freq * sqrt(number of spots) 
        high frequecy FMs will be calculated. If 'infer', the ratio_high_freq 
        will be set to 0. The default is 'infer'.
        A high can achieve better smothness. c should be setted to [0, 0.05].
    ratio_neighbors: float | 'infer', optional
        The ratio_neighbors will be used to determine the number of neighbors
        when contruct the KNN graph by spatial coordinates. Indeed, ratio_neig-
        hobrs * sqrt(number of spots) / 2 indicates the K. If 'infer', the para
        will be set to 1.0. The default is 'infer'.
    c: float, optional
        c balance the smoothness and difference with previous expresssion.
    spatial_info : list or tupple, optional
        The column names of spaital coordinates in adata.obs_names or key
        in adata.varm_keys() to obtain spatial information. The default
        is ['array_row', 'array_col'].
    normalize_lap : bool. optional
        Whether need to normalize the Laplcian matrix. The default is False.
        

    Returns
    -------
    count_matrix: DataFrame

    """
    import scipy.sparse as ss
    if ratio_low_freq == 'infer':
        if adata.shape[0] <= 800:
            num_low_frequency = min(15 * int(np.ceil(np.sqrt(adata.shape[0]))),
                                    adata.shape[0])
        elif adata.shape[0] <= 2000:
            num_low_frequency = 12 *  int(np.ceil(np.sqrt(adata.shape[0])))
        elif adata.shape[0] <= 10000:
            num_low_frequency = 10 *  int(np.ceil(np.sqrt(adata.shape[0])))
        else:
            num_low_frequency = 4 * int(np.ceil(np.sqrt(adata.shape[0])))
    else:
        num_low_frequency = int(np.ceil(np.sqrt(adata.shape[0]) * \
                                        ratio_low_freq))
    if ratio_high_freq == 'infer':
        num_high_frequency = 0 * int(np.ceil(np.sqrt(adata.shape[0])))
    else:
        num_high_frequency = int(np.ceil(np.sqrt(adata.shape[0]) * \
                                        ratio_high_freq))
        
    if ratio_neighbors == 'infer':
        if adata.shape[0] <= 500:
            num_neighbors = 4
        else:
            num_neighbors = int(np.ceil(np.sqrt(adata.shape[0]) / 2))
    else:
        num_neighbors = int(np.ceil(np.sqrt(adata.shape[0]) / 2 \
                                    * ratio_neighbors))
        
    adata.var_names_make_unique()
    sc.pp.filter_genes(adata, min_cells=1)
    # Get Laplacian matrix according to coordinates 
    lap_mtx = get_laplacian_mtx(adata, 
                                num_neighbors=num_neighbors,
                                spatial_key=spatial_info,
                                normalization=normalize_lap)
    
    # Fourier modes of low frequency
    eigvals_s, eigvecs_s = ss.linalg.eigsh(lap_mtx.astype(float), 
                                           k=num_low_frequency,
                                           which='SM')
    if num_high_frequency > 0:
        # Fourier modes of high frequency    
        eigvals_l, eigvecs_l = ss.linalg.eigsh(lap_mtx.astype(float),
                                               k=num_high_frequency, 
                                               which='LM')
        eigvals = np.concatenate((eigvals_s, eigvals_l))         # eigenvalues
        eigvecs = np.concatenate((eigvecs_s, eigvecs_l), axis=1) # eigenvectors
    else:
        eigvals = eigvals_s
        eigvecs = eigvecs_s
    
    # *********************** Graph Fourier Tranform **************************
    # Calculate GFT
    eigvecs_T = eigvecs.transpose()
    if not ss.issparse(adata.X):
        exp_mtx = adata.X
    else:
        exp_mtx = adata.X.toarray()
    frequency_array = np.matmul(eigvecs_T, exp_mtx)
    # low-pass filter
    filter_list = [1 / (1 + c * eigv) for eigv in eigvals]
    filter_array = np.matmul(np.diag(filter_list), frequency_array)
    filter_array = np.matmul(eigvecs, filter_array)
    filter_array[filter_array < 0] = 0
    
    # whether need to replace original count matrix
    if inplace and not ss.issparse(adata.X):
        adata.X = filter_array
    elif inplace:
        import scipy.sparse as ss
        adata.X = ss.csr.csr_matrix(filter_array)
        
    filter_array = pd.DataFrame(filter_array, 
                                index=adata.obs_names,
                                columns=adata.var_names)
    return filter_array

def select_num_fms(adata,
                   ratio_low_freq='infer', 
                   ratio_high_freq='infer', 
                   ratio_neighbors='infer',
                   spatial_info=['array_row', 'array_col'],
                   select_auto=True,
                   cutoff_fms = 0.05,
                   normalized_lap=True):
    """
    Select FMs automatically acoording to corresponding frequencies.

    Parameters
    ----------
    adata : AnnData
        adata.X is the normalized count matrix. Besides, the spatial coordinat-
        es of all spots could be found in adata.obs or adata.obsm.
    ratio_low_freq : float | "infer", optional
        The ratio_low_freq will be used to determine the number of the FMs with
        low frequencies. Indeed, the ratio_low_freq * sqrt(number of spots) low
        frequecy FMs will be calculated. If 'infer', the ratio_low_freq will be
        set to 1.0. The default is 'infer'.
    ratio_high_freq: float | 'infer', optional
        The ratio_high_freq will be used to determine the number of the FMs with
        high frequencies. Indeed, the ratio_high_freq * sqrt(number of spots) 
        high frequecy FMs will be calculated. If 'infer', the ratio_high_freq 
        will be set to 0. The default is 'infer'.
        A high can achieve better smothness. c should be setted to [0, 0.05].
    ratio_neighbors: float | 'infer', optional
        The ratio_neighbors will be used to determine the number of neighbors
        when contruct the KNN graph by spatial coordinates. Indeed, ratio_neig-
        hobrs * sqrt(number of spots) / 2 indicates the K. If 'infer', the para
        will be set to 1.0. The default is 'infer'.
    spatial_info : list or tupple, optional
        The column names of spaital coordinates in adata.obs_names or key
        in adata.varm_keys() to obtain spatial information. The default
        is ['array_row', 'array_col'].
    select_auto : bool, optional
        Determine the number of FMs automatically.. The default is True.
    cutoff_fms : float, optional
        Amount of information retained. The default is 0.95.
    normalized_lap : TYPE, optional
        DESCRIPTION. The default is True.

    Raises
    ------
    ValueError
        cutoff_fms should be in (0, 1]

    Returns
    -------
    float | DataFrame
        If select == True, return the number of FMs used, value * sqrt(N).
        Otherwise, return a dataframe contains infromation mount.

    """
    # Ensure parameters
    if ratio_low_freq == 'infer':
        num_low_frequency = int(np.ceil(np.sqrt(adata.shape[0])))
    else:
        num_low_frequency = int(np.ceil(np.sqrt(adata.shape[0]) * \
                                        ratio_low_freq))
    if ratio_high_freq == 'infer':
        num_high_frequency = int(np.ceil(np.sqrt(adata.shape[0])))
    else:
        num_high_frequency = int(np.ceil(np.sqrt(adata.shape[0]) * \
                                        ratio_high_freq))
        
    if ratio_neighbors == 'infer':
        num_neighbors = int(np.ceil(np.sqrt(adata.shape[0]) / 2))
    else:
        num_neighbors = int(np.ceil(np.sqrt(adata.shape[0]) / 2 \
                                    * ratio_neighbors))
    if adata.shape[0] <= 500:
        num_neighbors = 4
    
    # *************** Construct graph and corresponding matrixs ***************
    lap_mtx = get_laplacian_mtx(adata, num_neighbors=num_neighbors,
                                spatial_key=spatial_info,
                                normalization=normalized_lap)
    
    # Next, calculate the eigenvalues and eigenvectors of the Laplace matrix
    if num_high_frequency > 0:
        # Fourier bases of low frequency
        eigvals_s, eigvecs_s = ss.linalg.eigsh(lap_mtx.astype(float), 
                                               k=num_low_frequency,
                                               which='SM')  
        # Fourier bases of high frequency
        eigvals_l, eigvecs_l = ss.linalg.eigsh(lap_mtx.astype(float),
                                               k=num_high_frequency, 
                                               which='LM')       
        adata.uns['low_freq'] = eigvals_s
        adata.uns['high_freq'] = eigvals_l
        adata.uns['low_fms'] = eigvecs_s
        adata.uns['high_fms'] = eigvecs_l
    else:
        # Fourier bases of low frequency
        eigvals_s, eigvecs_s = ss.linalg.eigsh(lap_mtx.astype(float), 
                                               k=num_low_frequency,
                                               which='SM') 
        adata.uns['low_freq'] = eigvals_s
        adata.uns['low_fms'] = eigvecs_s
    
    # ********************** select fms ********************************
    eigvals_s = np.abs(eigvals_s)
    eigvals_s_power =  2 ** (-2 * eigvals_s)
    if select_auto:
        if cutoff_fms <= 0 or cutoff_fms > 1:
            raise ValueError("cutoff_fms should be in (0, 1]")
        frac_d = np.sum(eigvals_s_power)
        condition = 1
        for num in range(eigvals_s.size):
            if condition <= 1 - cutoff_fms:
                break
            condition -= eigvals_s_power[num] / frac_d
        
        return num / np.sqrt(adata.shape[1])
    else:
        plot_df = pd.DataFrame({'FMs': range(1, eigvals_s.size + 1),
                                'amount of information': eigvals_s_power})
        
        return plot_df
    
def rank_gene_smooth(adata, 
                     ratio_low_freq='infer', 
                     ratio_high_freq='infer', 
                     ratio_neighbors='infer',
                     spatial_info=['array_row', 'array_col'],
                     cal_pval=True,
                     normalize_lap=False,
                     exp_scale=True,
                     filter_peaks=False,
                     S=10,
                     spec_norm=True,
                     n_low_freq=None,
                     a=2):
    """
    Rank genes to find spatially variable genes by graph Fourier transform.

    Parameters
    ----------
    adata : AnnData
        adata.X is the normalized count matrix. Besides, the spatial coordinat-
        es could be found in adata.obs or adata.obsm.
    ratio_low_freq : float | "infer", optional
        The ratio_low_freq will be used to determine the number of the FMs with
        low frequencies. Indeed, the ratio_low_freq * sqrt(number of spots) low
        frequecy FMs will be calculated. If 'infer', the ratio_low_freq will be
        set to 1.0. The default is 'infer'.
    ratio_high_freq: float | 'infer', optional
        The ratio_high_freq will be used to determine the number of the FMs with
        high frequencies. Indeed, the ratio_high_freq * sqrt(number of spots) 
        high frequecy FMs will be calculated. If 'infer', the ratio_high_freq 
        will be set to 1.0. The default is 'infer'.
    ratio_neighbors: float | 'infer', optional
        The ratio_neighbors will be used to determine the number of neighbors
        when contruct the KNN graph by spatial coordinates. Indeed, ratio_neig-
        hobrs * sqrt(number of spots) / 2 indicates the K. If 'infer', the para
        will be set to 1.0. The default is 'infer'.
    spatial_name : list or tupple, optional
        The column names of spaital coordinates in adata.obs_names. The default
        is ['array_row', 'array_col'].
    spatial_key : string or None,
        If spatial_keq is provided, the spatial information could be found in
        adata.obsm['spaital_key'] (spot * coordinate). The default if none.
    cal_pval : bool, optional
        Whether need to calculate p val by mannwhitneyu. The default is False.
    normalize_lap : bool, optional
        Whether need to normalize laplacian matrix. The default is false.
    exp_scale : bool, optional
        If False, the expression matrix will not be scaled. The default is 
        True.
    filter_peaks: bool, optional
        For calculated vectors/signals in frequency/spectral domian, whether
        filter low peaks to stress the important peaks. The default is False.
    q : float, optional
        Quantile or sequence of quantiles to compute in frequency/spectral 
        domian. The default is 0.2.
    spec_norm : None or {'l1', 'l2', 'max'}, optional
        If None, the signal in spectral domian matrix will not be normzalized. 
        Otherwise, the given spec_norm method will be implemented. The default
        is l1. 
    Returns
    -------
    score_df : dataframe
        Return gene information.

    """
    # Ensure parameters
    if ratio_low_freq == 'infer':
        num_low_frequency = int(np.ceil(np.sqrt(adata.shape[0])))
    else:
        num_low_frequency = int(np.ceil(np.sqrt(adata.shape[0]) * \
                                        ratio_low_freq))
    if ratio_high_freq == 'infer':
        num_high_frequency = int(np.ceil(np.sqrt(adata.shape[0])))
    else:
        num_high_frequency = int(np.ceil(np.sqrt(adata.shape[0]) * \
                                        ratio_high_freq))
        
    if ratio_neighbors == 'infer':
        num_neighbors = int(np.ceil(np.sqrt(adata.shape[0]) / 2))
    else:
        num_neighbors = int(np.ceil(np.sqrt(adata.shape[0]) / 2 \
                                    * ratio_neighbors))
    if adata.shape[0] <= 500:
        num_neighbors = 4
            
    # Ensure gene index uniquely and all gene had expression  
    adata.var_names_make_unique()
    sc.pp.filter_genes(adata, min_cells=1)
    
    # *************** Construct graph and corresponding matrixs ***************
    if 'low_freq' not in adata.uns_keys():
        lap_mtx = get_laplacian_mtx(adata, num_neighbors=num_neighbors,
                                    spatial_key=spatial_info,
                                    normalization=normalize_lap)
        
        # Next, calculate the eigenvalues and eigenvectors of the Laplace matrix
        # Fourier bases of low frequency
        eigvals_s, eigvecs_s = ss.linalg.eigsh(lap_mtx.astype(float), 
                                               k=num_low_frequency,
                                               which='SM')  
        if n_low_freq != None:
            eigvals_s = eigvals_s[:n_low_freq]
            eigvecs_s = eigvecs_s[:n_low_freq, :]
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
    else:
        eigvals_s = adata.uns['low_freq']
        eigvecs_s = adata.uns['low_fms']
        if n_low_freq != None:
            eigvals_s = eigvals_s[:n_low_freq]
            eigvecs_s = eigvecs_s[:n_low_freq, :]
        if 'high_freq' in adata.uns_keys():
            eigvals_l = adata.uns['high_freq']
            eigvecs_l = adata.uns['high_fms']
            eigvals = np.concatenate((eigvals_s, eigvals_l))
            eigvecs = np.concatenate((eigvecs_s, eigvecs_l), axis=1)
        else:
            eigvals = eigvals_s
            eigvecs = eigvecs_s
        
    # ************************ Graph Fourier Tranform *************************
    # Calculate GFT
    eigvecs_T = eigvecs.transpose()
    if exp_scale:
        if type(adata.X) == np.ndarray:
            exp_mtx = preprocessing.scale(adata.X)
        else:
            exp_mtx = preprocessing.scale(adata.X.toarray())

    else:
        if type(adata.X) == np.ndarray:
            exp_mtx = adata.X
        else:
            exp_mtx = adata.X.toarray()
        
    frequency_array = np.matmul(eigvecs_T, exp_mtx)
    frequency_array = np.abs(frequency_array)
    
    # Filter noise peaks
    if filter_peaks == True:
        frequency_array_thres_low = np.quantile(frequency_array[:num_low_frequency, :],
                                            q=0.5, axis=0)
        frequency_array_thres_high = np.quantile(frequency_array[num_low_frequency:, :],
                                            q=0.5, axis=0)
        # frequency_array_thres = np.mean(frequency_array, axis=0)
        for j in range(frequency_array.shape[1]):
            # tmp = frequency_array[:num_low_frequency, j]
            frequency_array[:num_low_frequency, :][frequency_array[:num_low_frequency, j] <= \
                            frequency_array_thres_low[j], j]= 0
            frequency_array[num_low_frequency:, :][frequency_array[num_low_frequency:, j] <= \
                            frequency_array_thres_high[j], j]= 0

    if spec_norm != None:
        frequency_array = preprocessing.normalize(frequency_array, 
                                                  norm='l1',
                                                  axis=0)

    eigvals = np.abs(eigvals)
    eigvals_power = a ** ( -eigvals * 2) 
    score_list = np.matmul(eigvals_power, frequency_array)
    score_max = np.matmul(eigvals_power, (1 / len(eigvals)) * \
                          np.ones(len(eigvals)))       
    score_list = score_list / score_max
    print("Graph Fourier Transform finished!")

    # Rank genes according to smooth score
    adata.var["gft_score"] = score_list 
    score_df =adata.var["gft_score"]
    score_df =pd.DataFrame(score_df)
    score_df = score_df.sort_values(by="gft_score", ascending=False) 
    score_df.loc[:, "svg_rank"] = range(1, score_df.shape[0] + 1)
    adata.var["svg_rank"] = score_df.reindex(adata.var_names).loc[:, "svg_rank"]
    print("SVG ranking could be found in adata.obs['Rank']")
    
    # Cutoff of gft_score
    from kneed import KneeLocator
    magic = KneeLocator(score_df.svg_rank.values, 
                        score_df.gft_score.values,
                        direction='decreasing',
                        curve='convex',
                        S=S)
    score_df['cutoff_gft_score'] = False
    score_df['cutoff_gft_score'][:(magic.elbow + 1)] = True
    adata.var['cutoff_gft_score'] = score_df['cutoff_gft_score']
    print("The spatially variable genes judged by gft_score could be found " +
          "in adata.obs['cutoff_gft_score']")
    adata.varm['freq_domain_svg'] = frequency_array.transpose()
    print("Gene signals in frequency domain could be found in " + 
          "adata.varm['freq_domain']")
    adata.uns['frequencies_svg'] = eigvals
    adata.uns['fms_low'] = eigvecs_s
    adata.uns['fms_high'] =eigvecs_l
    
    # ****************** calculate pval ***************************
    if cal_pval == True:
        # pval_list = _significant_test_permutation(exp_mtx=exp_mtx.transpose(),
        #                             gene_score=score_list, eigvals=eigvals,
        #                             eigvecs_T=eigvecs_T)
        pval_list = _test_significant_freq(freq_array=adata.varm['freq_domain_svg'],
                                           cutoff=num_low_frequency)
        from statsmodels.stats.multitest import multipletests
        qval_list = multipletests(np.array(pval_list), method='fdr_by')[1]
        adata.var['pvalue'] = pval_list
        adata.var['qvalue'] = qval_list
        score_df = adata.var.loc[score_df.index, :].copy()

    return score_df

def calculate_frequcncy_domain(adata, 
                               ratio_low_freq='infer', 
                               ratio_high_freq='infer', 
                               ratio_neighbors='infer',
                               spatial_info=['array_row', 'array_col'],
                               return_freq_domain=True,
                               normalize_lap=False,
                               exp_norm=True,
                               filter_peaks=True,
                               spec_norm=True):
    """
    Obtain gene signals in frequency/spectral domain for given genes in 
    adata.var_names.

    Parameters
    ----------
    adata : AnnData
        adata.X is the normalized count matrix. Besides, the spatial coordinat-
        es could be found in adata.obs or adata.obsm.
    ratio_low_freq : float | "infer", optional
        The ratio_low_freq will be used to determine the number of the FMs with
        low frequencies. Indeed, the ratio_low_freq * sqrt(number of spots) low
        frequecy FMs will be calculated. If 'infer', the ratio_low_freq will be
        set to 1.0. The default is 'infer'.
    ratio_high_freq: float | 'infer', optional
        The ratio_high_freq will be used to determine the number of the FMs with
        high frequencies. Indeed, the ratio_high_freq * sqrt(number of spots) 
        high frequecy FMs will be calculated. If 'infer', the ratio_high_freq 
        will be set to 1.0. The default is 'infer'.
    ratio_neighbors: float | 'infer', optional
        The ratio_neighbors will be used to determine the number of neighbors
        when contruct the KNN graph by spatial coordinates. Indeed, ratio_neig-
        hobrs * sqrt(number of spots) / 2 indicates the K. If 'infer', the para
        will be set to 1.0. The default is 'infer'.
    spatial_info : list | tupple | str, optional
        The column names of spaital coordinates in adata.obs_keys() or 
        key in adata.obsm_keys. The default is ['array_row','array_col'].
    return_freq_domain : bool, optional
        Whether need to return gene signals in frequency domain. The default is 
        True.
    normalize_lap : bool, optional
        Whether need to normalize laplacian matrix. The default is false.
    exp_norm: bool, optional
        Whether need to normalize count matrix to z-score. The default is True.
    filter_peaks: bool, optional
        For calculated vectors/signals in frequency/spectral domian, whether
        filter low peaks to stress the important peaks. The default is False.
    spec_norm: bool, optional
        If False, the signal in spectral domian matrix will not be normzalized. 
        Otherwise, the given l1_norm method will be implemented. The default
        is True

    Returns
    -------
    DataFrame, the index indicates the gene and the columns indicates correspo-
    nding frequecies/smoothness. 

    """ 
    # Critical parameters
    # Ensure parameters
    if ratio_low_freq == 'infer':
        num_low_frequency = int(np.ceil(np.sqrt(adata.shape[0])))
    else:
        num_low_frequency = int(np.ceil(np.sqrt(adata.shape[0]) * \
                                        ratio_low_freq))
    if ratio_high_freq == 'infer':
        num_high_frequency = int(np.ceil(np.sqrt(adata.shape[0])))
    else:
        num_high_frequency = int(np.ceil(np.sqrt(adata.shape[0]) * \
                                        ratio_high_freq))
    if adata.shape[0] >= 10000:
        num_high_frequency = int(np.ceil(np.sqrt(adata.shape[0])))
        
    if ratio_neighbors == 'infer':
        num_neighbors = int(np.ceil(np.sqrt(adata.shape[0]) / 2))
    else:
        num_neighbors = int(np.ceil(np.sqrt(adata.shape[0]) / 2 \
                                    * ratio_neighbors))
    if adata.shape[0] <= 500:
        num_neighbors = 4
    
    # Ensure gene index uniquely and all gene had expression
    adata.var_names_make_unique()
    sc.pp.filter_genes(adata, min_cells=1)
    
    # ************** Construct graph and corresponding matrixs ****************
    lap_mtx = get_laplacian_mtx(adata, num_neighbors=num_neighbors,
                                spatial_key=spatial_info,
                                normalization=normalize_lap)
    
    # Calculate the eigenvalues and eigenvectors of the Laplace matrix
    np.random.seed(123)
    if num_high_frequency > 0: 
        # Fourier modes of low frequency
        eigvals_s, eigvecs_s = ss.linalg.eigsh(lap_mtx.astype(float), 
                                               k=num_low_frequency,
                                               which='SM')  
        # Fourier modes of high frequency    
        eigvals_l, eigvecs_l = ss.linalg.eigsh(lap_mtx.astype(float),
                                               k=num_high_frequency, 
                                               which='LM')
        eigvals = np.concatenate((eigvals_s, eigvals_l))         # eigenvalues
        eigvecs = np.concatenate((eigvecs_s, eigvecs_l), axis=1) # eigenvectors
    else:
        eigvals_s, eigvecs_s = ss.linalg.eigsh(lap_mtx.astype(float), 
                                               k=num_low_frequency,
                                               which='SM')
        eigvecs = eigvecs_s
        eigvals = eigvals_s
    
    # ************************Graph Fourier Tranform***************************
    # Calculate GFT
    eigvecs = eigvecs.transpose()
    if exp_norm == True:
        if not ss.issparse(adata.X):
            exp_mtx = preprocessing.scale(adata.X)
        else:
            exp_mtx = preprocessing.scale(adata.X.toarray())
    else:
        if not ss.issparse(adata.X):
            exp_mtx = adata.X
        else:
            exp_mtx = adata.X.toarray()
    frequency_array = np.matmul(eigvecs, exp_mtx)
    frequency_array = np.abs(frequency_array)
    # Filter noise peaks
    if filter_peaks == True:
        frequency_array_thres = np.mean(frequency_array, axis=0)
        for j in range(adata.shape[0]):
            frequency_array[frequency_array[:, j] <= \
                            frequency_array_thres[j], j] = 0
    # Spectral domian normalization
    if spec_norm == True:
        frequency_array = preprocessing.normalize(frequency_array, 
                                                  norm='l1',
                                                  axis=0)
    
    # ********************** Results of GFT ***********************************
    frequency_df = pd.DataFrame(frequency_array, columns=adata.var_names, 
                            index=['low_spec_' + str(low) \
                                  for low in range(1, num_low_frequency + 1)] \
                                + ['high_spec_' + str(high) \
                                 for high in range(1, num_high_frequency + 1)])
    adata.varm['freq_domain'] = frequency_df.transpose()
    print("Gene signals in frequency domain could be found in " + \
         " adata.varm['freq_domain']")
    adata.uns['frequencies'] = eigvals
    
    if return_freq_domain:
        return frequency_df

def find_tissue_module(adata, 
                       n_genes='infer', 
                       ratio_fms='infer',
                       ratio_neighbors='infer',
                       spatial_info=['array_row', 'array_col'],
                       n_neighbors=15,
                       resolution=1,
                       random_state=0,
                       **kargs):
    

    # Find tissue module by grouping Spatially variable genes with similar 
    # spatial patterns acoording to louvain algorithm.
    # Check conditions and determine parameters
    if 'svg_rank' not in adata.var.columns:
        assert KeyError("adata.var['svg_rank'] is not available. Please run\
                        SpaGFT.rank_gene_smooth(adata) firstly.")
    if ratio_fms == 'infer':
        if adata.shape[0] <= 500:
            low_ratio_freq = 4
        elif adata.shape[0] <= 10000:
            low_ratio_freq = 2
        else:
            low_ratio_freq = 1
    if n_genes == 'infer':
        n_genes = adata.var_names[adata.var['cutoff_gft_score']].size
    elif not isinstance(n_genes, int):
        raise ValueError("n_genes should be int.")
    gene_score = adata.var.sort_values(by='svg_rank')

    tmp_adata = adata.copy()
    tmp_adata.X = adata.raw[:, adata.var_names].X
    sc.pp.log1p(tmp_adata)
    calculate_frequcncy_domain(tmp_adata, 
                               ratio_low_freq=ratio_fms,
                               ratio_high_freq=0,
                               ratio_neighbors=ratio_neighbors,
                               spatial_info=spatial_info,
                               filter_peaks=False,
                               return_freq_domain=False)
    # Create new anndata to store freq domain information
    gft_adata = sc.AnnData(tmp_adata.varm['freq_domain'])
    sc.pp.neighbors(gft_adata, n_neighbors=n_neighbors, use_rep='X')
    sc.tl.umap(gft_adata)
    adata.varm['gft_umap'] = gft_adata.obsm['X_umap']
    # clustering
    gft_adata = gft_adata[gene_score.index[:n_genes], :]
    sc.pp.neighbors(gft_adata, n_neighbors=n_neighbors, use_rep='X')
    sc.tl.louvain(gft_adata, resolution=resolution, random_state=random_state,
                  **kargs)
    adata.var['tm_genes'] = 'None'
    adata.var.loc[gft_adata.obs_names, 'tm_genes'] = gft_adata.obs.louvain
    adata.var['tm_genes'] = pd.Categorical(adata.var['tm_genes'])
    # tm pseudo expression
    all_tms = gft_adata.obs.louvain.cat.categories
    tm_df = pd.DataFrame(0, index=adata.obs_names, columns='tm_' + all_tms)
    for tm in all_tms:
        pseudo_exp = np.ravel(tmp_adata[:,
                    gft_adata.obs.louvain[gft_adata.obs.louvain==tm].index].X.sum(axis=1))
        tm_df['tm_' + str(tm)] = pseudo_exp
    adata.obsm['tm_expression'] = tm_df
    tm_df = tm_df.copy()
    tm_df[tm_df < np.quantile(tm_df, q=0.85, axis=0)] = 0
    tm_df[tm_df >= np.quantile(tm_df, q=0.85, axis=0)] = 1
    adata.obsm['tm_region'] = tm_df
    # sub tm expression clustering
    tm_df = pd.DataFrame(index=adata.obs_names)
    for tm in all_tms:
        tm_gene_list = gft_adata.obs.louvain[gft_adata.obs.louvain==tm].index
        sub_gft_adata = gft_adata[tm_gene_list, :].copy()
        sc.pp.neighbors(sub_gft_adata, n_neighbors=n_neighbors, use_rep='X')
        sc.tl.louvain(sub_gft_adata, resolution=0.5, random_state=random_state,
                      **kargs)
        all_sub_tms = sub_gft_adata.obs.louvain.cat.categories
        for sub_tm in all_sub_tms:
            pseudo_exp = tmp_adata[:,
            sub_gft_adata.obs.louvain[sub_gft_adata.obs.louvain==sub_tm].index].X.sum(axis=1)
            tm_df['tm-' + str(tm) + "_subTm-" + str(sub_tm)] = pseudo_exp
    adata.obsm['subTm_expression'] = tm_df
    tm_df = tm_df.copy()
    tm_df[tm_df < np.quantile(tm_df, q=0.85, axis=0)] = 0
    tm_df[tm_df >= np.quantile(tm_df, q=0.85, axis=0)] = 1
    adata.obsm['subTm_region'] = tm_df

        
      
