# -*- coding: utf-8 -*-
# @Time : 2022/5/21 13:44
# @Author : Tory Deng
# @File : cluster.py
# @Software: PyCharm
from threading import Thread
from typing import Literal, Optional, Union, Callable
from multiprocessing.pool import ThreadPool

import pandas as pd
from scipy.spatial.distance import squareform
from functools import partial
import os

import anndata as ad
import leidenalg
import numpy as np
import scanpy as sc
from TracyWidom import TracyWidom
from loguru import logger
from sklearn.metrics.pairwise import paired_distances
from sklearn.mixture import GaussianMixture
from sklearn.neighbors import kneighbors_graph
from sklearn.preprocessing import minmax_scale
from sklearn.preprocessing import scale

from ..distance import compute_pairwise_distances


def do_clustering(
        adata: ad.AnnData,
        mode: Literal['one-way', 'two-way'],
        gene_clustering_method: Literal['gmm', 'leiden'],
        n_gene_clusters: Optional[int] = None,
        gene_distance: Union[str, Callable] = None,
        n_gene_neighbors: Optional[int] = 15,
        n_cell_clusters: Optional[int] = None,
        random_stat: Optional[int] = None
):
    """
    Only do gene clustering if `mode = 'one-way'`. Also simultaneously cluster cells and find highly confident cells
    if `mode = 'two-way'`.

    :param adata: The AnnData object
    :param mode: `one-way` only considers patterns in genes; `two-way` considers patterns both in cells and genes
    :param gene_clustering_method: The gene clustering method
    :param n_gene_clusters: The number of gene clusters
    :param gene_distance: The metric used to measure the distances between genes
    :param n_gene_neighbors: The number of gene neighbors
    :param n_cell_clusters: The number of cell clusters
    :param random_stat: The random seed. Default is none
    :return:
    """
    t1 = Thread(
        target=cluster_genes,
        args=(adata, mode, gene_clustering_method, n_gene_clusters, gene_distance, n_gene_neighbors, random_stat)
    )
    t1.start()
    if mode == 'two-way':
        t2 = Thread(
            target=find_high_confidence_cells,
            args=(adata, n_cell_clusters,),
            kwargs={'random_stat': random_stat}
        )
        t2.start()
    t1.join()
    if 't2' in locals().keys():
        t2.join()
        logger.debug(f"Size of cell clusters: \n{adata.obs['cluster'].value_counts()}")

    cluster_gene_counts = adata.var['cluster'].value_counts()
    logger.debug(f"Total number of gene clusters: {cluster_gene_counts.shape[0]}")
    logger.debug(f"Size of the top 5 largest gene clusters: \n{cluster_gene_counts.head(5)}")
    logger.debug(f"Size of the top 5 smallest gene clusters: \n{cluster_gene_counts.tail(5)}")


def cluster_genes(
        adata: ad.AnnData,
        mode: Literal['one-way', 'two-way'],
        method: Literal['gmm', 'leiden'],
        n_gene_clusters: Optional[int] = None,
        gene_distance: Union[str, Callable] = None,
        n_gene_neighbors: Optional[int] = 30,
        random_stat: Optional[int] = None,
):
    """
    Cluster genes using GMM or leiden clustering. For GMM clustering, the number of clusters is estimated when
    `n_gene_clusters` is None and `centrality` is computed when `mode` = 'one-way'. For leiden clustering, the SNN
    graph is calculated and the eigenvector centrality is also calculated when `mode` = 'one-way'.

    :param adata: The AnnData object
    :param mode: `one-way` only considers patterns in genes; `two-way` considers patterns both in cells and genes
    :param method: The gene clustering method
    :param n_gene_clusters: Number of gene clusters
    :param gene_distance: The metric used to measure the distances between genes
    :param n_gene_neighbors: Number of gene neighbors
    :param random_stat: The random seed. Default is none
    :return: None
    """
    logger.info(f"Start to cluster genes using {method} clustering...")

    if method == 'gmm':
        if n_gene_clusters is None:
            logger.info("You didn't specify the number of gene clusters. Automatically finding...")
            n_gene_clusters = estimate_K(adata.varm['pca'], pval=0.05, zero_replace=50)
        _check_gene_clustering_params(adata, method, n_gene_clusters, gene_distance)
        gmm = GaussianMixture(n_components=n_gene_clusters, init_params='k-means++', random_state=random_stat)
        adata.var['cluster'] = gmm.fit_predict(adata.varm['pca'])
        if mode == 'one-way':
            adata.var['centrality'] = _compute_distance2center(adata, gmm.means_)
    else:
        _check_gene_clustering_params(adata, method, n_gene_clusters, gene_distance)
        find_neighbors(adata, 'gene', gene_distance, n_gene_neighbors)
        if mode == 'one-way':
            adata.var['cluster'], adata.var['centrality'] = leiden(adata, 'gene', compute_centrality=True, seed=random_stat)
        else:
            adata.var['cluster'] = leiden(adata, 'gene', seed=random_stat)
    logger.info("Gene clustering finished!")


def find_high_confidence_cells(
        adata: ad.AnnData,
        n_cell_clusters: Optional[int] = None,
        high_prob: float = 0.95,
        high_freq: int = 8,
        min_cluster_size: int = 10,
        random_stat: Optional[int] = None,
):
    logger.info(f"Start to find highly confident cells...")
    # determine the number of cell clusters
    if n_cell_clusters is None:
        logger.info("You didn't specify the number of cell clusters. Automatically finding...")
        n_cell_clusters = estimate_K(adata.obsm['pca'])

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
            cluster_labels[~is_confident] = -1
            break
    adata.obs['cluster'], adata.obs['highly_confident'] = cluster_labels, is_confident
    logger.info(f"High-confidence cell detection finished!")


def estimate_K(X: np.ndarray, zero_replace: Union[int, str] = 'raise', pval: float = 0.001) -> int:
    """
    Estimate the optimal k for clustering algorithms. The function finds the eigenvalues of the sample covariance matrix.
    It will then return the number of significant eigenvalues according to the Tracy-Widom test.

    :param X: Processed input expression matrix of shape (n_cells, n_genes)
    :param zero_replace: 'raise' or int, default='raise'. Raise an error or replace the K when the estimated K is zero.
    :param pval: P-value
    :return: an estimated number of clusters k
    """
    tw = TracyWidom(beta=1)
    p, n = X.shape[0], X.shape[1]
    X_scaled = scale(X, axis=1)  # scale cells
    muTW = (np.sqrt(n - 1) + np.sqrt(p)) ** 2
    sigmaTW = (np.sqrt(n - 1) + np.sqrt(p)) * (1 / np.sqrt(n - 1) + 1 / np.sqrt(p)) ** (1 / 3)
    sigmaHatNaive = np.dot(X_scaled, X_scaled.T)
    bd = tw.cdfinv(1.0 - pval) * sigmaTW + muTW

    evals = np.linalg.eigvalsh(sigmaHatNaive)
    k = int((evals > bd).sum())
    if k == 0:
        if zero_replace == 'raise':
            raise RuntimeError("The estimated number of clusters equals 0!")
        else:
            logger.warning(f"Replacing the estimated number of clusters (0) with {zero_replace}...")
            k = zero_replace
    else:
        logger.info(f"Estimated number of clusters: {k}")
    return k


def find_neighbors(
        adata: ad.AnnData,
        on: Literal['cell', 'gene'],
        dis_metric: Union[str, Callable] = 'euclidean',
        n_neighbors: int = 15,
):
    attr = 'obs' if on == 'cell' else 'var'
    attrp, attrm = attr + 'p', attr + 'm'

    logger.debug(f"Compute {on} distances...")
    getattr(adata, attrp)['distances'] = compute_pairwise_distances(
        getattr(adata, attrm)['pca'], metric=dis_metric, n_jobs=-1
    )

    logger.debug(f"Finding KNNs on {on}...")
    getattr(adata, attrp)['connectivities'] = kneighbors_graph(
        getattr(adata, attrp)['distances'], n_neighbors, metric='precomputed', n_jobs=-1
    ).toarray()


def leiden(data: Union[ad.AnnData, np.ndarray],
           on: Optional[Literal['cell', 'gene']] = None,
           resolution: float = 1.0,
           compute_centrality: bool = False,
           seed: Optional[int] = None
           ):
    if isinstance(data, ad.AnnData):
        if on is None:
            raise TypeError("`on` can not be None when `data` is `ad.AnnData` object!")
        attr = 'obs' if on == 'cell' else 'var'
        attrp, attrm = attr + 'p', attr + 'm'
        if 'connectivities' not in getattr(data, attrp):
            raise ValueError(f"Can not find `connectivities` in {attrp}. Please calculate neighbors first.")
        adjacency = getattr(data, attrp)['connectivities']
    elif isinstance(data, np.ndarray):
        if data.ndim != 2:
            raise ValueError(f"Dimension of data must be 2, but got {data.ndim}!")
        if data.shape[0] != data.shape[1]:
            raise ValueError(f"`data` should be a square matrix!")
        adjacency = data
    else:
        raise TypeError(f"`data` can only be `ad.AnnData` or `np.ndarray`!")

    directed = False if np.allclose(adjacency, adjacency.T) else True
    logger.debug(f"Creating {'directed' if directed else 'undirected'} graph...")
    G = sc._utils.get_igraph_from_adjacency(adjacency, directed=directed)
    logger.debug("Leiden clustering starts...")
    partition = leidenalg.find_partition(G, partition_type=leidenalg.RBConfigurationVertexPartition, weights=G.es['weight'],
                                         n_iterations=-1, resolution_parameter=resolution, seed=seed)
    cluster_labels = np.array(partition.membership)
    if compute_centrality:
        logger.debug("Computing centrality...")
        centralities = np.zeros_like(cluster_labels)
        for cluster in np.unique(cluster_labels):
            cluster_mask = cluster_labels == cluster
            G_sub = G.subgraph(np.argwhere(cluster_mask).squeeze())
            centralities[cluster_mask] = G_sub.evcent()
        return cluster_labels, centralities
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
                dense_co_membership.append(1.0)
            else:
                dense_co_membership.append(0.0)
    return np.array(dense_co_membership)


def _check_gene_clustering_params(
        adata: ad.AnnData,
        method: Literal['gmm', 'leiden'],
        n_gene_clusters: Optional[int],
        gene_distance: Union[str, Callable]
):
    if method == 'gmm':
        if isinstance(n_gene_clusters, int):
            if n_gene_clusters < 2 or n_gene_clusters > adata.n_vars:
                raise ValueError(f"Argument `n_gene_clusters` must be between 2 and the number of genes in GMM.")
        else:
            raise TypeError(f"Argument `n_gene_clusters` must be integer or None, not {type(n_gene_clusters)} in GMM.")
        if gene_distance is not None:
            raise ValueError("Don't need to specify distance metric in GMM clustering!")
    elif method == 'leiden':
        if n_gene_clusters is not None:
            raise ValueError(f"Can not set `n_gene_clusters` in leiden clustering.")
        if gene_distance is None:
            raise ValueError("Must specify distance metric to create graph in leiden clustering!")
    else:
        raise NotImplementedError(f"`{method}` has not been implemented yet.")


def _compute_distance2center(adata: ad.AnnData, means_: np.ndarray) -> np.ndarray:
    """
    Compute the distance of each gene to its cluster mean, normalized in cluster.

    :param adata: The AnnData object
    :param means_: The cluster means
    :return: distances of all genes
    """
    all_distances = paired_distances(adata.varm['pca'], means_[adata.var['cluster']])
    for gene_cluster in np.unique(adata.var['cluster']):
        gene_cluster_mask = adata.var['cluster'] == gene_cluster
        all_distances[gene_cluster_mask] = minmax_scale(all_distances[gene_cluster_mask])
    return all_distances
