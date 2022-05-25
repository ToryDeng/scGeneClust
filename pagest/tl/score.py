# -*- coding: utf-8 -*-
# @Time : 2022/5/22 23:43
# @Author : Tory Deng
# @File : score.py
# @Software: PyCharm
import os
from multiprocessing.pool import ThreadPool
from typing import Literal

import anndata as ad
import numpy as np
from loguru import logger
from scipy.spatial.distance import squareform
from scipy.stats import kruskal
from scipy.stats import spearmanr
from sklearn.feature_selection import f_classif
from sklearn.metrics import silhouette_samples


def score_gene_cluster(adata: ad.AnnData, method: Literal['silhouette', 'spearman']):
    """
    Compute scores of all gene clusters.

    :param adata: The AnnData object
    :param method: 'silhouette': mean silhouette scores; 'spearman': mean absolute spearman correlation
    :return: None
    """
    logger.info(f"Start to score gene clusters using {method}...")
    if method == 'silhouette':
        if 'distance' in adata.varp:
            sample_scores = silhouette_samples(adata.varp['distance'], adata.var['cluster'], metric='precomputed')
        else:
            sample_scores = silhouette_samples(adata.varm['svd'], adata.var['cluster'], metric='euclidean')
        for cluster in adata.var['cluster'].unique():
            cluster_mask = adata.var['cluster'] == cluster
            adata.var.loc[cluster_mask, 'cluster_score'] = sample_scores[cluster_mask].mean()
    elif method == 'spearman':
        corr = spearmanr(adata.varm['svd'], axis=1)[0]
        for cluster in adata.var['cluster'].unique():
            cluster_mask = adata.var['cluster'] == cluster
            cluster_corr = corr[np.ix_(cluster_mask, cluster_mask)]
            adata.var.loc[cluster_mask, 'cluster_score'] = np.mean(np.abs(squareform(cluster_corr, checks=False)))
    else:
        raise NotImplementedError(f"`{method}` has not been implemented yet.")
    logger.info(f"Scores of gene clusters have been saved in adata.var['cluster_score']")


def score_discriminative_gene(adata: ad.AnnData, stat: Literal['kw', 'f'], use_rep: str = 'log-normalized'):
    """
    Compute the Kruskal–Wallis H test statistic or F statistic of each gene after filtering low-quality gene clusters.
    This function is only called when the `mode` is 'one-way'.

    :param adata: The AnnData object
    :param stat: The kind of statistic used to represent differences of gene expression between cell groups
    :param use_rep: If use_rep is None, the statistic is calculated on adata.X; Otherwise calculated on adata.layers[use_rep]
    :return: None
    """
    X, y = adata.X if use_rep is None else adata.layers[use_rep], adata.obs['cluster']
    if stat == 'f':
        adata.var['stat'], _ = f_classif(X, y)
    elif stat == 'kw':
        # split expression of cell clusters
        cclusters = [X[y == grp, :] for grp in y.unique()]
        pool = ThreadPool(processes=os.cpu_count() - 1)
        result = pool.starmap(kruskal, [[expr[:, i] for expr in cclusters] for i in range(adata.n_vars)], chunksize=100)
        adata.var['stat'] = [res[0] for res in result]  # stat, pval
    else:
        raise NotImplementedError(f"{stat} has not been implemented yet!")
