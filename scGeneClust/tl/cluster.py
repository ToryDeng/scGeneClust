# -*- coding: utf-8 -*-
# @Time : 2022/6/3 17:01
# @Author : Tory Deng
# @File : cluster.py
# @Software: PyCharm
from threading import Thread
from typing import Literal, Optional
from multiprocessing.pool import ThreadPool

import pandas as pd
from scipy.spatial.distance import squareform
from functools import partial
import os

import anndata as ad
import numpy as np
from loguru import logger
from sklearn.mixture import GaussianMixture
from ._utils import _compute_1w_gene_score, _find_neighbors, leiden, _compute_cell_co_membership


def do_clustering(
        adata: ad.AnnData,
        mode: Literal['one-way', 'two-way'],
        n_gene_clusters: Optional[int] = None,
        n_gene_neighbors: Optional[int] = 30,
        n_cell_clusters: Optional[int] = None,
        random_stat: Optional[int] = None,
        **kwargs
):
    t1 = Thread(
        target=cluster_genes,
        args=(adata, mode, n_gene_clusters, n_gene_neighbors, random_stat)
    )
    t1.start()

    if mode == 'two-way':
        high_prob = kwargs.get('high_prob', 0.95)
        high_freq = kwargs.get('high_freq', 8)
        min_cluster_size = kwargs.get('min_cluster_size', 10)

        t2 = Thread(
            target=find_high_confidence_cells,
            args=(adata, n_cell_clusters, high_prob, high_freq, min_cluster_size, random_stat),
        )
        t2.start()

    t1.join()
    if 't2' in locals().keys():
        t2.join()


def cluster_genes(
        adata: ad.AnnData,
        mode: Literal['one-way', 'two-way'],
        n_gene_clusters: Optional[int],
        n_gene_neighbors: Optional[int],
        random_stat: Optional[int],
):
    if mode == 'one-way':
        gmm = GaussianMixture(n_components=n_gene_clusters, init_params='k-means++', random_state=random_stat)
        adata.var['cluster'] = gmm.fit_predict(adata.varm['pca'])
        adata.var['score'] = _compute_1w_gene_score(adata, gmm.means_)
    else:
        _find_neighbors(adata, 'gene', 'spearman', n_gene_neighbors)
        adata.var['cluster'] = leiden(adata, 'gene', resolution=5, seed=random_stat)


def find_high_confidence_cells(
        adata: ad.AnnData,
        n_cell_clusters: int,
        high_prob: float = 0.95,
        high_freq: int = 8,
        min_cluster_size: int = 10,
        random_stat: Optional[int] = None,
):
    logger.info(f"Start to find high-confidence cells...")

    # compute the frequency matrix
    pool = ThreadPool(processes=os.cpu_count() - 1)
    top_comps = [adata.obsm['pca'][:, :i] for i in range(1, 11)]
    partial_compute = partial(_compute_cell_co_membership, n_clusters=n_cell_clusters, seed=random_stat, p=high_prob)
    results = pool.map(partial_compute, top_comps)
    frequency_matrix = squareform(np.sum(results, axis=0))

    for freq_th in range(high_freq, 0, -1):
        cut_matrix = np.where(frequency_matrix < freq_th, 0, frequency_matrix)
        cluster_labels = leiden(cut_matrix, seed=random_stat)
        cluster_counts = pd.Series(cluster_labels).value_counts(ascending=False)
        cut_k = max(n_cell_clusters, np.argwhere((cluster_counts < min_cluster_size).values).squeeze()[0])
        valid_clusters = cluster_counts.index[:cut_k]
        is_confident = np.isin(cluster_labels, valid_clusters)
        if is_confident.sum() / is_confident.shape[0] > 0.05:
            logger.debug(f"Final frequency cutoff: {freq_th}")
            cluster_labels[~is_confident] = -1
            break
    adata.obs['cluster'], adata.obs['highly_confident'] = cluster_labels, is_confident
    logger.info(f"High-confidence cell detection finished!")
