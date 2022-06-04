# -*- coding: utf-8 -*-
# @Time : 2022/5/25 19:41
# @Author : Tory Deng
# @File : distance.py
# @Software: PyCharm
from typing import Union, Optional, Literal

import numpy as np
import anndata as ad
from scipy.sparse import csr_matrix
from sklearn.metrics.pairwise import pairwise_distances, _VALID_METRICS
from scipy.spatial.distance import squareform
from scipy.stats import spearmanr


def compute_pairwise_distances(
        adata: ad.AnnData,
        on: Literal['cell', 'gene'],
        metric: str = 'rho_p',
        n_jobs: Optional[int] = None,
        **kwargs
) -> Union[np.ndarray, csr_matrix]:
    """
    Compute pairwise distances between cells or genes.

    :param adata: The AnnData object
    :param on: compute distance matrix on cells or genes
    :param metric: The distance metric
    :param n_jobs: The number of jobs to use for the computation
    :param kwargs: Other keyword arguments
    :return: ndarray or csr_matrix, distance matrix
    """
    if metric in _VALID_METRICS and metric != 'correlation':
        X = get_array(adata, use_rep='pca', use_emb=True, on=on)
        return pairwise_distances(X=X, metric=metric, n_jobs=n_jobs, **kwargs)
    elif metric == 'spearman':
        X = get_array(adata, use_rep='pca', use_emb=True, on=on)
        corr, _ = spearmanr(X, axis=1)
        np.fill_diagonal(corr, 1)
        return 1 - np.abs(corr)
    elif metric == 'pearson':
        X = get_array(adata, use_rep='pca', use_emb=True, on=on)
        corr = np.corrcoef(X)
        np.fill_diagonal(corr, 1)
        return 1 - np.abs(corr)
    elif metric == 'rho_p':
        X = get_array(adata, use_rep=f'X_{on}_log', use_emb=False, on=on)
        data = X.T + 1 / (X.shape[0] ** 2)
        var_log = np.var(np.log(data), axis=0)
        var_mtx = [var_log[i] + var_log[j] for i in range(data.shape[1] - 1) for j in range(i + 1, (data.shape[1]))]
        var_mtx = squareform(var_mtx) + np.diag(2 * var_log)
        return 2 * np.cov(np.log(data.T), ddof=0) / var_mtx
    elif metric == 'phi_s':
        X = get_array(adata, use_rep=f'X_{on}_log', use_emb=False, on=on)
        data = X.T + 1 / X.shape[0] ** 2
        var_log = np.var(np.log(data), axis=0)
        phi_s = [np.var(np.log(data[:, i]) - np.log(data[:, j])) / var_log[i] for i in range(data.shape[1]) for j in
                 range(data.shape[1])]
        phi_s = -1.0 * np.array(phi_s).reshape(data.shape[1], data.shape[1])
        return phi_s
    else:
        raise NotImplementedError(f"`{metric}` has not been implemented yet.")


def get_array(
        adata: ad.AnnData,
        use_rep: Optional[str],
        use_emb: bool,
        on: Literal['cell', 'gene'],
):
    # prepare data matrix, in which rows correspond to obs, and cols correspond to vars
    if use_emb:
        X = getattr(adata, 'obsm' if on == 'cell' else 'varm')[use_rep]
    else:
        X = adata.layers[use_rep] if on == 'cell' else adata.layers[use_rep].T
    return X







