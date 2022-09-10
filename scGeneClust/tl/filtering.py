# -*- coding: utf-8 -*-
# @Time : 2022/6/3 17:55
# @Author : Tory Deng
# @File : filtering.py
# @Software: PyCharm
from typing import Literal, Optional

import anndata as ad
import numpy as np
from loguru import logger
from sklearn.ensemble import IsolationForest
from sklearn.feature_selection import mutual_info_classif


def handle_single_gene_cluster(adata: ad.AnnData, version: Literal['fast', 'ps'], random_stat: Optional[int]):
    gene_cluster_counts = adata.var['cluster'].value_counts()
    # summary of gene gene_clustering_graph
    logger.debug(f"Total number of gene clusters: {gene_cluster_counts.shape[0]}")
    logger.debug(f"Total number of single gene clusters: {(gene_cluster_counts == 1).sum()}")
    bins = [0.001, 1, 5, 1e1, 5e1, 1e2, 5e2, 1e3, 5e3, 1e4, 5e4]
    cluster_size_counts = gene_cluster_counts.value_counts(bins=bins, sort=False)
    logger.debug(f"Bins of the gene cluster size: \n{cluster_size_counts[cluster_size_counts > 0]}")
    is_single_cluster = adata.var['cluster'].isin(gene_cluster_counts[gene_cluster_counts == 1].index)

    if is_single_cluster.sum() > 0:
        if version == 'fast':
            single_cluster_genes = adata.var_names[is_single_cluster]
            deviances = compute_deviance(adata.X[:, is_single_cluster])
            is_outlier = IsolationForest(random_state=random_stat).fit_predict(deviances.reshape(-1, 1)) == -1
            keep_genes = np.logical_or(~is_single_cluster, adata.var_names.isin(single_cluster_genes[is_outlier]))
            logger.debug(f"Removing {adata.n_vars - keep_genes.sum()} single gene clusters...")
            adata._inplace_subset_var(keep_genes)
        else:
            adata.var['representative'] = False
            grouped = adata.var[~is_single_cluster].groupby(by='cluster')['score']
            for i in range(len(grouped.nlargest(1))):
                adata.var.loc[grouped.nlargest(1).index[i][1], 'representative'] = True
            logger.debug(f"Number of representative genes in non-single cluster: {adata.var['representative'].sum()}")
            thres = adata.var.loc[adata.var['representative'], 'score'].min()
            adata.var.loc[is_single_cluster, 'representative'] = adata.var.loc[is_single_cluster, 'score'] > thres
    else:
        logger.debug("Not found any single gene cluster!")


def filter_constant_genes(adata: ad.AnnData):
    is_constant_genes = np.all(adata.X == adata.X[0, :], axis=0)
    logger.debug(f"Removing {is_constant_genes.sum()} constant genes...")
    adata._inplace_subset_var(~is_constant_genes)


def filter_low_confidence_cells(adata: ad.AnnData):
    logger.debug(f"Size of cell clusters: \n{adata.obs['cluster'].value_counts()}")
    logger.debug(f"Removing {adata.n_obs - adata.obs['highly_confident'].sum()} low-confidence cells...")
    adata._inplace_subset_obs(adata.obs['highly_confident'])


def compute_deviance(X: np.ndarray):
    """Compute deviance score for each single gene cluster"""
    pi = X.sum(0) / X.sum()
    n = X.sum(1)[:, np.newaxis]
    with np.errstate(all='ignore'):
        left_half = np.nansum(X * np.log(X / n.dot(pi[np.newaxis, :])), axis=0)
        right_half = np.nansum((n - X) * (np.log((n - X) / (n.dot((1 - pi)[np.newaxis, :])))), axis=0)
    return 2 * (left_half + right_half)


def filter_irrelevant_gene(adata: ad.AnnData, top_pct: int, random_stat: int):
    """Find relevant genes according to mutual information with cluster labels of highly confident cells"""
    logger.info(f"Start to find relevant genes...")
    relevance = mutual_info_classif(adata.layers['X_gene_log'], adata.obs.cluster,
                                    discrete_features=False, random_state=random_stat)
    adata.var['score'] = relevance
    rlv_th = np.percentile(relevance, 100 - top_pct if top_pct is not None else 80)
    logger.debug(f"Relevance threshold: {rlv_th}")
    adata._inplace_subset_var(relevance > rlv_th)
    logger.debug(f"Candidate relevant genes: {adata.shape[1]}")
    logger.info(f"Relevant gene detection finished!")
