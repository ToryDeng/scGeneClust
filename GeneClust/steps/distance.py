# -*- coding: utf-8 -*-
# @Time : 2022/4/8 21:01
# @Author : Tory Deng
# @File : distance.py
# @Software: PyCharm
"""
Metrics used to calculate the distance between genes
"""
from typing import Literal

import anndata as ad
import numpy as np
import pandas as pd
from scipy.spatial.distance import pdist, squareform
from sklearn.metrics import pairwise_distances
from sklearn.feature_selection import mutual_info_regression
import compositional as comp


def compute_gene_distance(
        adata: ad.AnnData,
        metric: Literal['pearson', 'spearman', 'kendall', 'bayesian', 'euclidean', 'mahalanobis', 'rho_p', 'phi_s']
):
    """
    Compute the similarity matrix for genes and store it in adata.varp.

    :param adata: The anndata object
    :param metric: The similarity metric
    :return: The similarity matrix, in which each entry is the similarity between two genes
    """
    if metric in ('pearson', 'spearman', 'kendall'):  # cols represent genes after transpose
        adata.varp[metric] = 1 - adata.varm[adata.uns['dr_method']].T.corr(metric).abs()
    elif metric == 'bayesian':
        adata.varp[metric] = 1 - bayes_corr(pd.DataFrame(adata.raw.X.T, index=adata.raw.var_names))
    elif metric == 'euclidean':
        adata.varp[metric] = euclidean_dis(adata.varm[adata.uns['dr_method']])
    elif metric == 'mahalanobis':
        adata.varp[metric] = mahalanobis_dis(adata.varm[adata.uns['dr_method']])
    elif metric == 'rho_p':
        adata.varp[metric] = rho_p(adata)
    elif metric == 'phi_s':
        adata.varp[metric] = phi_s(adata)
    else:
        raise NotImplementedError(f"Metric {metric} has not been implemented!")
    adata.uns['distance'] = metric
    print(f"distances between genes have been computed using {metric}.")


def bayes_corr(data: pd.DataFrame) -> pd.DataFrame:
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
    # calculate mutual_info between different genes(pairwise)
    simi_matrix = squareform(
        pdist(data, lambda u, v: mutual_info_regression(np.array(u).reshape(-1, 1), v, random_state=40))
    )
    simi_matrix = pd.DataFrame(simi_matrix, columns=data.index, index=data.index)
    # calculate mutual_info between a gene and itself
    s = [float(mutual_info_regression(np.array(data[i]).reshape(-1, 1), data[i])) for i in range(len(data))]
    simi_matrix = simi_matrix + np.diag(s)

    return simi_matrix


def euclidean_dis(data: pd.DataFrame):
    """
    similarity measure using euclidean distance
    :param data: row: gene
    """
    dist = pd.DataFrame(squareform(pdist(data, 'euclidean')), columns=data.index, index=data.index)
    return dist


def mahalanobis_dis(data: pd.DataFrame):
    """
    similarity measure using mahalanobis distance
    :param data: rows: genes
    """
    dist = pd.DataFrame(squareform(pdist(data, 'mahalanobis')), columns=data.index, index=data.index)
    return dist


def rho_p(adata: ad.AnnData):
    """
    similarity measure using rho_p
    :param data: rows: genes
    Citation:
    * https://github.com/tpq/propr
    * https://www.nature.com/articles/s41592-019-0372-4
    equation: rho_p = 2cov(logx, logy)/(var(logx)+var(logy))
    """
    data = adata.raw.X + 1/(adata.raw.X.shape[1]**2)
    # propr_matrix = pairwise_distances(data.T, metric=lambda u, v: 2*np.cov(np.log(u),np.log(v), ddof=0)[0][1]/(np.var(np.log(u)) + np.var(np.log(v))), n_jobs=-1)
    # return pd.DataFrame(propr_matrix, columns=adata.raw.var_names, index=adata.raw.var_names)
    var_log = np.var(np.log(data), axis=0)
    var_mtx = [var_log[i]+var_log[j] for i in range(data.shape[1] - 1) for j in range(i + 1, (data.shape[1]))]
    var_mtx = squareform(var_mtx) + np.diag(2*var_log)

    return pd.DataFrame(2*np.cov(np.log(data.T), ddof=0) / var_mtx, columns=adata.raw.var_names, index=adata.raw.var_names)


def phi_s(adata: ad.AnnData):
    """
    similarity measure using phi_s
    :param data: rows: genes
    Citation:
    * https://github.com/tpq/propr
    * https://www.nature.com/articles/s41592-019-0372-4
    equation: phi_s(log(x), log(y)) = var(log(x/y)) /  var(log(x))
    """
    data = adata.raw.X + 1/adata.raw.X.shape[1]**2
    var_log = np.var(np.log(data), axis=0)
    phi_s = [np.var(np.log(data[:,i])-np.log(data[:,j])) / var_log[i] for i in range(data.shape[1]) for j in range(data.shape[1])]
    phi_s = np.array(phi_s).reshape(data.shape[1], data.shape[1]) 
    return pd.DataFrame(phi_s, columns=adata.raw.var_names, index=adata.raw.var_names)
    # propr_matrix = pairwise_distances(data.T, metric=lambda u, v: np.var(np.log(u) - np.log(v)) / np.var(np.log(u)), n_jobs=-1)
    # return pd.DataFrame(propr_matrix, columns=adata.raw.var_names, index=adata.raw.var_names)







