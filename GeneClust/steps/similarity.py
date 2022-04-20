# -*- coding: utf-8 -*-
# @Time : 2022/4/8 21:01
# @Author : Tory Deng
# @File : similarity.py
# @Software: PyCharm
"""
Metrics used to calculate the similarity between genes
"""
import pandas as pd
import numpy as np
import anndata as ad
from sklearn.feature_selection import mutual_info_regression
from scipy.spatial.distance import pdist, squareform
from typing import Literal


def compute_gene_similarity(
        adata: ad.AnnData,
        dr_method: str,
        metric: Literal[
            'pearson', 'spearman', 'kendall', 'bayes_corr', 'mutual_info', 'euclidean_dis', 'mahalanobis_dis']
):
    """
    Compute the similarity matrix for genes and store it in adata.varp.

    :param adata: The anndata object
    :param dr_method: The dimension reduction method
    :param metric: The similarity metric
    :return: The similarity matrix, in which each entry is the similarity between two genes
    """
    if metric in ('pearson', 'spearman', 'kendall'):
        adata.varp[metric] = adata.varm[dr_method].T.corr(metric).abs()  # the cols represent genes after transpose
    elif metric == 'bayes_corr':
        adata.varp[metric] = bayes_corr(pd.DataFrame(adata.raw.X.T, index = adata.raw.var_names))
    elif metric == 'mutual_info':
        adata.varp[metric] = mutual_info(adata.varm[dr_method])
    elif metric == 'euclidean_dis':
        adata.varp[metric] = euclidean_dis(adata.varm[dr_method])
    elif metric == 'mahalanobis_dis':
        adata.varp[metric] = mahalanobis_dis(adata.varm[dr_method])
    else:
        raise NotImplementedError(f"Metric {metric} has not been implemented!")


def bayes_corr(data: pd.DataFrame):
    """
    similarity measure using Bayesian correlation 
    :param data: raw data (row:gene) 
    """
    nrowsX = data.shape[0]
    ncolsX = data.shape[1]
    alpha0 = [1 / nrowsX] * ncolsX
    beta0 = [1 - x for x in alpha0]
    cs = data.sum(axis=0).tolist()
    alphas = np.array(data) + alpha0
    betas = np.array(beta0 * nrowsX).reshape(nrowsX, -1) + np.array(cs * nrowsX).reshape(nrowsX, -1) - np.array(data)
    alphasPLUSbetas = alphas + betas
    alp_alpPLUSbeta = alphas / alphasPLUSbetas
    Psi = alp_alpPLUSbeta - np.array([x / ncolsX for x in alp_alpPLUSbeta.sum(axis=1)] * ncolsX).reshape(-1, nrowsX).T
    var_vec = ((((alphas * betas) / ((alphasPLUSbetas ** 2) * (alphasPLUSbetas + 1))).sum(axis=1) 
                + (Psi ** 2).sum(axis=1)) / ncolsX).reshape(nrowsX, 1)
    cov_mtrx = np.dot(Psi, Psi.T) / ncolsX
    Bcorrvals = cov_mtrx / np.sqrt(np.dot(var_vec, var_vec.T)) 
    Bcorrvals[np.diag_indices_from(Bcorrvals)] = 1

    return pd.DataFrame(Bcorrvals, columns=data.index, index=data.index)

def mutual_info(data: pd.DataFrame):
    """
    similarity measure using mutual information
    :param data: row: gene
    """
    #calculate mutual_info between different genes(pairwise)
    simi_matrix = squareform(pdist(data, lambda u, v: mutual_info_regression(np.array(u).reshape(-1, 1), v, random_state = 40)))
    simi_matrix = pd.DataFrame(simi_matrix, columns = data.index, index = data.index)
    #calculate mutual_info between a gene and itself
    s = [float(mutual_info_regression(np.array(data[i]).reshape(-1,1),data[i])) for i in range(len(data))]
    simi_matrix = simi_matrix + np.diag(s)

    return simi_matrix
    # In each iteration, calculate the mi between i_th gene and all genes
    #data_array shape(ngene, ncell)
    #data_array = data.values
    # for i in range(len(data)):
        # simi_matrx[i] = mutual_info_regression(data_array.T, data.iloc[i,], random_state = 40)
        #Or try this 
        # simi_matrix.iloc[i:, i] = mutual_info_regression(data_array[i:].T, data.iloc[i,], random_state = 40)
        # t = simi_matrix.T.copy()
        # np.fill_diagonal(t.values, 0)
        # simi_matrix = simi_matrix + t


def euclidean_dis(data: pd.DataFrame):
    """
    similarity measure using euclidean distance
    :param data: row: gene
    """
    dist = pd.DataFrame(squareform(pdist(data, 'euclidean')), columns = data.index, index = data.index)
    return dist


def mahalanobis_dis(data: pd.DataFrame):
    """
    similarity measure using mahalanobis distance
    :param data: rows: genes
    """
    dist = pd.DataFrame(squareform(pdist(data, 'mahalanobis')), columns = data.index, index = data.index)
    return dist