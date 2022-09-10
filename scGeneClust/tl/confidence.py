# -*- coding: utf-8 -*-
# @Time : 2022/6/3 17:01
# @Author : Tory Deng
# @File : confidence.py
# @Software: PyCharm
from functools import partial
from multiprocessing import cpu_count
from multiprocessing.pool import ThreadPool
from typing import Literal, Optional, Union

import anndata as ad
import leidenalg
import numpy as np
import pandas as pd
import scanpy as sc
from loguru import logger
from scipy.spatial.distance import squareform
from sklearn.mixture import GaussianMixture


def leiden(
        data: Union[ad.AnnData, np.ndarray],
        on: Literal['cell', 'gene'] = None,
        resolution: float = 1.0,
        seed: Optional[int] = None
):
    if isinstance(data, ad.AnnData):
        if on is None:
            raise TypeError("`on` can not be None when `data` is `ad.AnnData` object!")
        attr = 'obs' if on == 'cell' else 'var'
        attrp, attrm = attr + 'p', attr + 'm'
        if 'connectivities' not in getattr(data, attrp):
            raise ValueError(f"Can not find `connectivities` in {attrp}. Please calculate {on} neighbors first.")
        adjacency = getattr(data, attrp)['connectivities']
    elif isinstance(data, np.ndarray):
        if data.ndim != 2:
            raise ValueError(f"Dimension of data must be 2, but got {data.ndim} dim(s) array!")
        if data.shape[0] != data.shape[1]:
            raise ValueError(f"`data` should be a square matrix!")
        adjacency = data
        on = 'given matrix'
    else:
        raise TypeError(f"`data` can only be `ad.AnnData` or `np.ndarray`!")

    directed = False if np.allclose(adjacency, adjacency.T) else True
    logger.debug(f"Creating {'directed' if directed else 'undirected'} graph on {on}...")
    G = sc._utils.get_igraph_from_adjacency(adjacency, directed=directed)
    logger.debug("Leiden gene_clustering_graph starts...")
    partition = leidenalg.find_partition(G, partition_type=leidenalg.RBConfigurationVertexPartition,
                                         weights=G.es['weight'],
                                         n_iterations=-1, resolution_parameter=resolution, seed=seed)
    cluster_labels = np.array(partition.membership)
    logger.debug("Leiden gene_clustering_graph finished!")
    return cluster_labels


def _compute_cell_co_membership(
        X: np.ndarray,
        n_clusters: int,
        seed: Optional[int] = None,
        p: float = 0.95
) -> np.ndarray:
    X = X.reshape(-1, 1) if X.ndim == 1 else X
    gmm = GaussianMixture(n_components=n_clusters, init_params='k-means++', random_state=seed)
    gmm.fit(X)
    cluster_labels = gmm.predict(X)
    dense_co_membership = []
    cell_proba = gmm.predict_proba(X)
    max_proba = cell_proba.max(1)
    for cell_i in range(X.shape[0] - 1):
        for cell_j in range(cell_i + 1, X.shape[0]):
            if cluster_labels[cell_i] == cluster_labels[cell_j] and max_proba[cell_i] > p and max_proba[cell_j] > p:
                dense_co_membership.append(1.)
            else:
                dense_co_membership.append(0.)
    return np.array(dense_co_membership)


def find_high_confidence_cells(
        adata: ad.AnnData,
        n_cell_clusters: int,
        random_stat: Optional[int] = None,
):
    logger.info(f"Finding high-confidence cells...")

    # compute the frequency matrix
    pool = ThreadPool(processes=cpu_count() - 1)
    top_comps = [adata.obsm['pca'][:, :i] for i in range(2, 12)]
    partial_compute = partial(_compute_cell_co_membership, n_clusters=n_cell_clusters, seed=random_stat, p=0.95)
    results = pool.map(partial_compute, top_comps)
    frequency_matrix = squareform(np.sum(results, axis=0))

    for freq_th in range(8, 0, -1):
        cut_matrix = np.where(frequency_matrix < freq_th, 0, frequency_matrix)
        cluster_labels = leiden(cut_matrix, seed=random_stat, on='cell')
        cluster_counts = pd.Series(cluster_labels).value_counts(ascending=False)
        cut_k = min(n_cell_clusters, np.argwhere((cluster_counts > 10).values).squeeze()[0])
        valid_clusters = cluster_counts.index[:cut_k]
        is_confident = np.isin(cluster_labels, valid_clusters)
        if is_confident.sum() / is_confident.shape[0] > 0.05:
            logger.debug(f"Final frequency cutoff: {freq_th}")
            cluster_labels[~is_confident] = -1
            break
    adata.obs['cluster'], adata.obs['highly_confident'] = cluster_labels, is_confident
    logger.info(f"High-confidence cell detection finished!")
