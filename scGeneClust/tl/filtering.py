# -*- coding: utf-8 -*-
# @Time : 2022/6/3 17:55
# @Author : Tory Deng
# @File : filtering.py
# @Software: PyCharm
from typing import Literal,Optional

import anndata as ad
import numpy as np
from loguru import logger
from sklearn.ensemble import IsolationForest
from sklearn.feature_selection import mutual_info_classif


def handle_single_gene_cluster(adata: ad.AnnData, mode: Literal['fast', 'hc'], random_stat: Optional[int]):
    gene_cluster_counts = adata.var['cluster'].value_counts()
    # summary of gene clustering
    logger.debug(f"Total number of gene clusters: {gene_cluster_counts.shape[0]}")
    logger.debug(f"Total number of single gene clusters: {(gene_cluster_counts == 1).sum()}")
    bins = [0.001, 1, 5, 1e1, 5e1, 1e2, 5e2, 1e3, 5e3, 1e4, 5e4]
    cluster_size_counts = gene_cluster_counts.value_counts(bins=bins, sort=False)
    logger.debug(f"Bins of the gene cluster size: \n{cluster_size_counts[cluster_size_counts > 0]}")

    
    # filter_adata inliers in single gene clusters
    if mode == 'fast':
        is_single_cluster = adata.var['cluster'].isin(gene_cluster_counts[gene_cluster_counts == 1].index)
        if is_single_cluster.sum() > 0:
            single_cluster_genes = adata.var_names[is_single_cluster]
            deviances = compute_deviance(adata.X[:, is_single_cluster])
            is_outlier = IsolationForest(random_state=random_stat).fit_predict(deviances.reshape(-1, 1)) == -1
            keep_genes = np.logical_or(~is_single_cluster, adata.var_names.isin(single_cluster_genes[is_outlier]))
            logger.debug(f"Removing {adata.n_vars - keep_genes.sum()} single gene clusters...")
            adata._inplace_subset_var(keep_genes)
    else:
        thres = min(adata.var.loc[adata.var['representative'] == True, 'score'])
        for i in adata.var[adata.var['cluster'] == -1].index:
            if adata.var['score'][i] > thres:
                    adata.var['representative'][i] = True
        adata._inplace_subset_var(adata.var['representative'])

def filter_constant_genes(adata: ad.AnnData):
    is_constant_genes = np.all(adata.X == adata.X[0, :], axis=0)
    logger.debug(f"Removing {is_constant_genes.sum()} constant genes...")
    adata._inplace_subset_var(~is_constant_genes)


def filter_low_confidence_cells(adata: ad.AnnData):
    logger.debug(f"Size of cell clusters: \n{adata.obs['cluster'].value_counts()}")
    logger.debug(f"Removing {adata.n_obs - adata.obs['highly_confident'].sum()} low-confidence cells...")
    adata._inplace_subset_obs(adata.obs['highly_confident'])


def compute_deviance(X: np.ndarray):
    pi = X.sum(0) / X.sum()
    n = X.sum(1)[:, np.newaxis]
    with np.errstate(all='ignore'):
        left_half = np.nansum(X * np.log(X / n.dot(pi[np.newaxis, :])), axis=0)
        right_half = np.nansum((n - X) * (np.log((n - X) / (n.dot((1 - pi)[np.newaxis, :])))), axis=0)
    return 2 * (left_half + right_half)


def filter_unrelevant_gene(
        adata: ad.AnnData,
        rlv_threshold: float = 0.01,
        random_stat: int = 42):

    #find relevant gene according to mutual information with cluster label based on high confident cells
    logger.info(f"Start to find relevant genes...")
    relevance = mutual_info_classif(adata.layers['X_gene_log'],adata.obs.cluster,discrete_features=True,random_state=random_stat)
    adata.var['score'] = relevance
    adata._inplace_subset_var(relevance > rlv_threshold)
    logger.info(f"Relevant gene detection finished!")